import streamlit as st
import google.generativeai as genai
from datetime import date, timedelta
import pandas as pd
import gspread
from gspread_dataframe import get_as_dataframe
from Bio import Entrez
import time

# --- HẰNG SỐ GIỚI HẠN ---
HISTORY_LIMIT = 40
COLLECTION_LIMIT = 5
DRUGS_PER_COLLECTION_LIMIT = 7
PRESCRIPTION_LIMIT_PER_DAY = 5

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
PROMPT_PRESCRIPTION = load_prompt("prompt_prescription.txt")

@st.cache_resource
def get_lookup_model():
    model_name = st.secrets.get("models", {}).get("lookup", "gemini-1.5-flash-latest")
    return genai.GenerativeModel(model_name)

@st.cache_resource
def get_pro_model():
    model_name = st.secrets.get("models", {}).get("pro", "gemini-pro")
    return genai.GenerativeModel(model_name)

@st.cache_resource
def get_prescription_model():
    model_name = st.secrets.get("models", {}).get("prescription", "gemini-1.5-pro-latest")
    return genai.GenerativeModel(model_name)


# --- HÀM XỬ LÝ GOOGLE SHEETS & PRO ---
@st.cache_data(ttl=600)
def get_access_codes_df():
    try:
        credentials_dict = dict(st.secrets.connections.gsheets.credentials)
        gspread_client = gspread.service_account_from_dict(credentials_dict)
        
        spreadsheet_name = st.secrets.connections.gsheets.spreadsheet
        spreadsheet = gspread_client.open(spreadsheet_name)
        worksheet = spreadsheet.sheet1
        return get_as_dataframe(worksheet)
    except Exception as e:
        # Trả lại thông báo lỗi thân thiện
        st.error(f"Lỗi kết nối tới Google Sheets.")
        return pd.DataFrame()

def verify_code(db, user_info, user_code):
    if not user_code:
        return False, "Vui lòng nhập mã truy cập."
    user_id = user_info['localId']
    token = user_info['idToken']
    codes_df = get_access_codes_df()
    if codes_df.empty:
        return False, "Không thể tải dữ liệu mã truy cập."
    codes_df.dropna(subset=['code'], inplace=True)
    codes_df['code'] = codes_df['code'].astype(str)
    if user_code not in codes_df['code'].values:
        return False, "Mã không hợp lệ hoặc không tìm thấy."
    claimed_codes = db.child("claimed_pro_codes").get(token=token).val() or {}
    if user_code in claimed_codes:
        if claimed_codes[user_code] != user_id:
            return False, "Mã này đã được sử dụng bởi một tài khoản khác."
        else:
            st.session_state.pro_access = True
            return True, "Bạn đã kích hoạt mã này trước đó. Chào mừng trở lại!"
    try:
        db.child("claimed_pro_codes").child(user_code).set(user_id, token=token)
        db.child("user_data").child(user_id).child("is_pro").set(True, token=token)
        st.session_state.pro_access = True
        return True, "Kích hoạt PRO thành công! Cảm ơn bạn đã ủng hộ dự án."
    except Exception as e:
        return False, f"Đã xảy ra lỗi trong quá trình kích hoạt: {e}"

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

# --- HÀM TRA CỨU THUỐC ---
def get_drug_info_from_api(drug_name, is_pro_user=False):
    identifier_model = get_lookup_model()
    prompt_nhan_dien_final = PROMPT_NHAN_DIEN.format(drug_name=drug_name)
    response_nhan_dien = identifier_model.generate_content(prompt_nhan_dien_final)
    response_text = response_nhan_dien.text
    try:
        hoat_chat_goc = response_text.split("Output:")[1].strip()
    except IndexError:
        hoat_chat_goc = response_text.strip()
    if hoat_chat_goc == "INVALID" or not hoat_chat_goc:
        return f"❌ Lỗi: '{drug_name}' không được nhận dạng.", None

    analysis_model = get_lookup_model()
    analysis_prompt = PROMPT_REGULAR
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
        
    return final_response, hoat_chat_goc

# --- CÁC HÀM TƯƠNG TÁC FIREBASE ---
def load_user_data(db, user_info):
    try:
        user_id = user_info['localId']
        token = user_info['idToken']
        data = db.child("user_data").child(user_id).get(token=token).val()
        if not data:
            return [], {}, False
        history = data.get("history", [])
        collections = data.get("collections", {})
        is_pro = data.get("is_pro", False)
        if is_pro:
            st.session_state.pro_access = True
        return history, collections
    except Exception:
        return [], {}

def load_user_result(db, user_info, drug_name):
    try:
        user_id = user_info['localId']
        token = user_info['idToken']
        result = db.child("user_data").child(user_id).child("results_cache").child(drug_name).get(token=token).val()
        return result
    except Exception:
        return None

def save_new_result(db, user_info, drug_name, result_text):
    try:
        user_id = user_info['localId']
        token = user_info['idToken']
        is_pro = db.child("user_data").child(user_id).child("is_pro").get(token=token).val() or False
        db.child("user_data").child(user_id).child("results_cache").child(drug_name).set(result_text, token=token)
        history = db.child("user_data").child(user_id).child("history").get(token=token).val() or []
        if drug_name in history:
            history.remove(drug_name)
        history.insert(0, drug_name)
        if not is_pro and len(history) > HISTORY_LIMIT:
            drug_to_delete = history.pop()
            db.child("user_data").child(user_id).child("results_cache").child(drug_to_delete).remove(token=token)
        db.child("user_data").child(user_id).child("history").set(history, token=token)
        return history
    except Exception as e:
        st.error(f"Lỗi khi lưu kết quả tra cứu: {e}")
        return None

def create_new_collection(db, user_info, collection_name):
    if not collection_name or collection_name.isspace():
        return False, "Tên bộ sưu tập không được để trống."
    try:
        user_id = user_info['localId']
        token = user_info['idToken']
        is_pro = db.child("user_data").child(user_id).child("is_pro").get(token=token).val() or False
        collections = db.child("user_data").child(user_id).child("collections").get(token=token).val() or {}
        if not is_pro and len(collections) >= COLLECTION_LIMIT and collection_name not in collections:
            return False, f"Đã đạt giới hạn {COLLECTION_LIMIT} bộ sưu tập."
        if collection_name in collections:
            return False, f"Bộ sưu tập '{collection_name}' đã tồn tại."
        collections[collection_name] = True
        db.child("user_data").child(user_id).child("collections").set(collections, token=token)
        return True, f"Đã tạo thành công bộ sưu tập '{collection_name}'."
    except Exception as e:
        return False, f"Đã xảy ra lỗi không xác định: {e}"

def add_drug_to_collection(db, user_info, collection_name, drug_name):
    try:
        user_id = user_info['localId']
        token = user_info['idToken']
        is_pro = db.child("user_data").child(user_id).child("is_pro").get(token=token).val() or False
        collections = db.child("user_data").child(user_id).child("collections").get(token=token).val() or {}
        if collection_name not in collections:
            return f"Lỗi: Không tìm thấy bộ sưu tập '{collection_name}'."
        drug_list_or_placeholder = collections.get(collection_name)
        if drug_list_or_placeholder is True:
            drug_list = []
        elif isinstance(drug_list_or_placeholder, list):
            drug_list = drug_list_or_placeholder
        else:
            return "Lỗi: Cấu trúc dữ liệu của bộ sưu tập không hợp lệ."
        if drug_name in drug_list:
            return f"'{drug_name}' đã có trong bộ sưu tập này."
        if not is_pro and len(drug_list) >= DRUGS_PER_COLLECTION_LIMIT:
            return f"Bộ sưu tập '{collection_name}' đã đầy (tối đa {DRUGS_PER_COLLECTION_LIMIT} thuốc)."
        drug_list.append(drug_name)
        collections[collection_name] = drug_list
        db.child("user_data").child(user_id).child("collections").set(collections, token=token)
        return f"Đã thêm '{drug_name}' vào '{collection_name}'."
    except Exception as e:
        return f"Đã xảy ra lỗi không xác định: {e}"

def delete_from_history(db, user_info, drug_name):
    try:
        user_id = user_info['localId']
        token = user_info['idToken']
        history = db.child("user_data").child(user_id).child("history").get(token=token).val() or []
        if drug_name in history:
            history.remove(drug_name)
            db.child("user_data").child(user_id).child("history").set(history, token=token)
        return True, f"Đã xóa '{drug_name}' khỏi lịch sử."
    except Exception as e:
        return False, f"Lỗi khi xóa khỏi lịch sử: {e}"

def delete_from_collection(db, user_info, collection_name, drug_name):
    try:
        user_id = user_info['localId']
        token = user_info['idToken']
        collections = db.child("user_data").child(user_id).child("collections").get(token=token).val() or {}
        if collection_name in collections and isinstance(collections[collection_name], list):
            drug_list = collections[collection_name]
            if drug_name in drug_list:
                drug_list.remove(drug_name)
                if not drug_list:
                    collections[collection_name] = True
                else:
                    collections[collection_name] = drug_list
                db.child("user_data").child(user_id).child("collections").set(collections, token=token)
                return True, f"Đã xóa '{drug_name}' khỏi '{collection_name}'."
        return False, "Không tìm thấy thuốc hoặc bộ sưu tập."
    except Exception as e:
        return False, f"Lỗi khi xóa khỏi bộ sưu tập: {e}"

def delete_collection(db, user_info, collection_name):
    try:
        user_id = user_info['localId']
        token = user_info['idToken']
        db.child("user_data").child(user_id).child("collections").child(collection_name).remove(token=token)
        return True, f"Đã xóa bộ sưu tập '{collection_name}'."
    except Exception as e:
        return False, f"Lỗi khi xóa bộ sưu tập: {e}"

# --- HÀM CHO PHÂN TÍCH ĐƠN THUỐC ---
def get_prescription_analysis(db, user_info, patient_context, prescription_text):
    try:
        user_id = user_info['localId']
        token = user_info['idToken']
        is_pro = st.session_state.get("pro_access", False)
        
        if not is_pro:
            today_str = date.today().isoformat()
            usage_ref = db.child("user_data").child(user_id).child("usage_counters").child("prescription_analysis")
            usage_data = usage_ref.get(token=token).val() or {"count": 0, "last_updated": ""}
            
            if usage_data.get("last_updated") != today_str:
                usage_data = {"count": 0, "last_updated": today_str}
            
            if usage_data.get("count", 0) >= PRESCRIPTION_LIMIT_PER_DAY:
                return f"❌ Bạn đã hết {PRESCRIPTION_LIMIT_PER_DAY} lượt phân tích miễn phí trong ngày hôm nay."

        model = get_prescription_model()
        prompt = PROMPT_PRESCRIPTION.format(patient_context=patient_context, prescription_text=prescription_text)
        response = model.generate_content(prompt)
        
        if not is_pro:
            usage_data["count"] += 1
            usage_ref.set(usage_data, token=token)
            
        return response.text
        
    except Exception as e:
        return f" Lỗi trong quá trình phân tích: {e}"
