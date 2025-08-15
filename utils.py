import streamlit as st
import google.generativeai as genai
from datetime import date, timedelta
import pandas as pd
import gspread
from gspread_dataframe import get_as_dataframe
from Bio import Entrez
import time

# --- CÁC HÀM PROMPT VÀ KHỞI TẠO MODEL ---
def load_prompt(file_path):
    try:
        with open(file_path, "r", encoding="utf-8") as f: return f.read()
    except FileNotFoundError:
        st.error(f"LỖI: Không tìm thấy file prompt tại '{file_path}'.")
        st.stop()

PROMPT_NHAN_DIEN = load_prompt("prompt_nhandien.txt")
PROMPT_REGULAR = load_prompt("prompt_regular.txt")
PROMPT_PRO = load_prompt("prompt_pro.txt")
PROMPT_SUMMARY = load_prompt("prompt_summary.txt")

@st.cache_resource
def get_regular_model():
    model_name = st.secrets.get("models", {}).get("regular", "gemini-1.5-flash-latest")
    return genai.GenerativeModel(model_name)

@st.cache_resource
def get_pro_model():
    model_name = st.secrets.get("models", {}).get("pro", "gemini-pro")
    return genai.GenerativeModel(model_name)

# --- CÁC HÀM XỬ LÝ GOOGLE SHEETS ---
@st.cache_data(ttl=600)
def get_access_codes_df():
    try:
        credentials = st.secrets.connections.gsheets.credentials
        gspread_client = gspread.service_account_from_dict(credentials)
        spreadsheet = gspread_client.open(st.secrets.connections.gsheets.spreadsheet)
        worksheet = spreadsheet.sheet1
        return get_as_dataframe(worksheet)
    except Exception as e:
        st.error(f"Lỗi kết nối tới Google Sheets.")
        return pd.DataFrame()

def verify_code(user_code):
    if not user_code: return False, "Vui lòng nhập mã truy cập."
    codes_df = get_access_codes_df()
    if codes_df.empty: return False, "Không thể tải dữ liệu mã truy cập."
    codes_df.dropna(subset=['code'], inplace=True)
    codes_df['code'] = codes_df['code'].astype(str)
    matched_code_series = codes_df[codes_df['code'] == user_code]
    if matched_code_series.empty: return False, "Mã không hợp lệ hoặc không tìm thấy."
    code_info = matched_code_series.iloc[0]
    code_type = code_info['type']
    if code_type == 'permanent':
        st.session_state.pro_access = True
        return True, f"Xác thực thành công! Chào mừng {code_info.get('owner', 'Pro User')}."
    if code_type == 'temporary':
        try:
            created_date = pd.to_datetime(code_info['created_at']).date()
            today = date.today()
            days_passed = (today - created_date).days
            if 0 <= days_passed <= 7:
                st.session_state.pro_access = True
                return True, f"Xác thực thành công! Mã còn hiệu lực {7 - days_passed} ngày."
            else: return False, "Mã tạm thời đã hết hạn."
        except Exception: return False, "Lỗi định dạng ngày tháng trong Google Sheet."
    return False, "Loại mã không xác định."

# --- CÁC HÀM XỬ LÝ DỮ LIỆU & API ---
@st.cache_data(ttl=3600)
def search_pubmed(drug_name):
    Entrez.email = "duocdien.ai.project@example.com"
    api_key = st.secrets.get("api_keys", {}).get("pubmed")
    if api_key: Entrez.api_key = api_key
    today = date.today()
    two_years_ago = today - timedelta(days=730)
    date_filter = f'AND ("{two_years_ago.strftime("%Y/%m/%d")}"[Date - Publication] : "{today.strftime("%Y/%m/%d")}"[Date - Publication])'
    search_term = f'"{drug_name}"[Title/Abstract] AND ("clinical trial"[Publication Type] OR "systematic review"[Publication Type]) {date_filter}'
    try:
        handle = Entrez.esearch(db="pubmed", term=search_term, retmax="5", sort="relevance")
        record = Entrez.read(handle)
        handle.close()
        id_list = record["IdList"]
        if not id_list: return "Không tìm thấy bài báo phù hợp nào trong 2 năm gần đây trên PubMed."
        handle = Entrez.efetch(db="pubmed", id=id_list, rettype="medline", retmode="text")
        records_text = handle.read()
        handle.close()
        context = ""
        articles = records_text.strip().split("\n\n")
        for article_text in articles:
            title = next((line[6:] for line in article_text.split('\n') if line.startswith("TI  - ")), "N/A")
            abstract = next((line[6:] for line in article_text.split('\n') if line.startswith("AB  - ")), "N/A")
            journal = next((line[6:] for line in article_text.split('\n') if line.startswith("JT  - ")), "N/A")
            pub_date = next((line[6:] for line in article_text.split('\n') if line.startswith("DP  - ")), "N/A")
            pmid = next((line[6:] for line in article_text.split('\n') if line.startswith("PMID- ")), "N/A")
            context += f"- Tiêu đề: {title}\n- Tạp chí: {journal}\n- Năm: {pub_date[:4]}\n- Tóm tắt: {abstract}\n- PMID: {pmid.strip()}\n\n"
            time.sleep(0.1)
        return context
    except Exception as e:
        return f"Đã xảy ra lỗi khi truy vấn API của PubMed: {e}"

@st.cache_data(ttl="6h")
def get_drug_info(drug_name, is_pro_user=False):
    identifier_model = get_regular_model()
    prompt_nhan_dien_final = PROMPT_NHAN_DIEN.format(drug_name=drug_name)
    response_nhan_dien = identifier_model.generate_content(prompt_nhan_dien_final)
    response_text = response_nhan_dien.text
    try:
        hoat_chat_goc = response_text.split("Output:")[1].strip()
    except IndexError:
        hoat_chat_goc = response_text.strip()
    if hoat_chat_goc == "INVALID" or not hoat_chat_goc:
        return f"❌ Lỗi: '{drug_name}' không được nhận dạng."
    
    analysis_model = get_pro_model() if is_pro_user else get_regular_model()
    analysis_prompt = PROMPT_PRO if is_pro_user else PROMPT_REGULAR
    
    generation_config = {"max_output_tokens": 8192, "temperature": 0.6}
    full_prompt = f"{analysis_prompt}\n\nHãy tra cứu và trình bày thông tin cho thuốc sau đây: **{hoat_chat_goc}**"
    
    response_phan_tich = analysis_model.generate_content(full_prompt, generation_config=generation_config)
    base_response_text = response_phan_tich.text
    final_response = f"✅ Hoạt chất đã nhận diện: **{hoat_chat_goc}**\n\n---\n\n{base_response_text}"

    if is_pro_user:
        section_11_content = "\n\n---\n\n**11. Phân tích các Nghiên cứu Lâm sàng nổi bật (trong 2 năm gần đây):**\n"
        try:
            with st.spinner("Người dùng Pro: Đang truy vấn API của PubMed..."):
                search_context = search_pubmed(hoat_chat_goc)
                summary_prompt_final = PROMPT_SUMMARY.format(drug_name=hoat_chat_goc, search_results=search_context)
                summary_model = get_pro_model()
                summary_response = summary_model.generate_content(summary_prompt_final, generation_config=generation_config)
                section_11_content += summary_response.text
        except Exception as e:
            st.warning(f"Lỗi khi xử lý thông tin từ PubMed: {e}")
            section_11_content += "Đã xảy ra lỗi khi cố gắng tóm tắt dữ liệu từ PubMed."
        final_response += section_11_content
        
    return final_response

# --- CÁC HÀM TƯƠNG TÁC FIREBASE (LỊCH SỬ) ---
def load_user_history(db, user_info):
    """Tải lịch sử tra cứu của người dùng từ Firebase."""
    try:
        user_id = user_info['localId']
        token = user_info['idToken']
        history = db.child("user_data").child(user_id).child("history").get(token=token).val()
        return history if history else []
    except Exception as e:
        # st.error("Lỗi khi tải lịch sử tra cứu.")
        # st.exception(e)
        return []

def save_drug_to_history(db, user_info, drug_name):
    """Lưu một thuốc vào lịch sử tra cứu của người dùng trên Firebase."""
    try:
        user_id = user_info['localId']
        token = user_info['idToken']
        current_history = load_user_history(db, user_info)
        
        if drug_name not in current_history:
            current_history.insert(0, drug_name)
            if len(current_history) > 20:
                current_history = current_history[:20]
            db.child("user_data").child(user_id).child("history").set(current_history, token=token)
    except Exception as e:
        st.warning("Lỗi khi lưu lịch sử tra cứu.")

# --- CÁC HÀM TƯƠNG TÁC FIREBASE (BỘ SƯU TẬP) ---
def load_user_collections(db, user_info):
    """Tải toàn bộ bộ sưu tập của người dùng từ Firebase."""
    try:
        user_id = user_info['localId']
        token = user_info['idToken']
        collections = db.child("user_data").child(user_id).child("collections").get(token=token).val()
        return collections if collections else {}
    except Exception as e:
        # st.error("Lỗi khi tải các bộ sưu tập.")
        return {}

def add_drug_to_collection(db, user_info, collection_name, drug_name):
    """Thêm một thuốc vào một bộ sưu tập cụ thể."""
    try:
        user_id = user_info['localId']
        token = user_info['idToken']
        drug_list = db.child("user_data").child(user_id).child("collections").child(collection_name).get(token=token).val()
        if drug_list is None:
            drug_list = []
        if drug_name not in drug_list:
            drug_list.append(drug_name)
            db.child("user_data").child(user_id).child("collections").child(collection_name).set(drug_list, token=token)
            return True
        return False
    except Exception as e:
        st.warning(f"Lỗi khi thêm thuốc vào bộ sưu tập '{collection_name}'.")
        return False

def create_new_collection(db, user_info, collection_name):
    """Tạo một bộ sưu tập mới (rỗng)."""
    if not collection_name or collection_name.isspace():
        return False, "Tên bộ sưu tập không được để trống."
    try:
        user_id = user_info['localId']
        token = user_info['idToken']
        existing_collections = load_user_collections(db, user_info)
        if collection_name in existing_collections:
            return False, f"Bộ sưu tập '{collection_name}' đã tồn tại."
        db.child("user_data").child(user_id).child("collections").child(collection_name).set([], token=token)
        return True, f"Đã tạo thành công bộ sưu tập '{collection_name}'."
    except Exception as e:
        st.warning(f"Lỗi khi tạo bộ sưu tập '{collection_name}'.")
        return False, "Đã xảy ra lỗi không xác định."
