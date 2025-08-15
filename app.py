import streamlit as st
import google.generativeai as genai
import google.generativeai.types as genai_types
import google.api_core.exceptions as ga_ex
import datetime
from datetime import date, timedelta
import pandas as pd
import gspread
from gspread_dataframe import get_as_dataframe
from gspread.exceptions import SpreadsheetNotFound
from Bio import Entrez
import time
# --- CÁC THƯ VIỆN MỚI CHO TÍNH NĂNG ĐĂNG NHẬP ---
import yaml
from yaml.loader import SafeLoader
import streamlit_authenticator as stauth

# --- KIỂM TRA TRẠNG THÁI BẢO TRÌ ---
is_maintenance = st.secrets.get("maintenance_mode", False)
if is_maintenance:
    st.set_page_config(page_title="Bảo trì", page_icon="🛠️")
    st.title("🛠️ Dược Điển AI đang được bảo trì")
    message = st.secrets.get("maintenance_message", "Ứng dụng đang được cập nhật. Vui lòng quay lại sau.")
    st.info(message)
    st.stop()

# --- 1. KHỞI TẠO TRẠNG THÁI PHIÊN ---
if 'history' not in st.session_state: st.session_state.history = []
if 'pro_access' not in st.session_state: st.session_state.pro_access = False

# --- 2. CẤU HÌNH VÀ TẢI PROMPTS ---
def load_prompt(file_path):
    try:
        with open(file_path, "r", encoding="utf-8") as f: return f.read()
    except FileNotFoundError:
        st.error(f"LỖI: Không tìm thấy file prompt tại '{file_path}'.")
        st.stop()
try:
    genai.configure(api_key=st.secrets["GOOGLE_API_KEY"])
except (FileNotFoundError, KeyError):
    st.error("LỖI: Vui lòng cấu hình GOOGLE_API_KEY trong secrets.toml.")
    st.stop()
PROMPT_NHAN_DIEN = load_prompt("prompt_nhandien.txt")
PROMPT_REGULAR = load_prompt("prompt_regular.txt")
PROMPT_PRO = load_prompt("prompt_pro.txt")
PROMPT_SUMMARY = load_prompt("prompt_summary.txt")

# --- 3. CÁC HÀM XỬ LÝ (GIỮ NGUYÊN) ---

# --- HÀM XỬ LÝ MÃ TRUY CẬP ---
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

# --- HÀM XỬ LÝ DƯỢC ĐIỂN ---
@st.cache_resource
def get_regular_model():
    model_name = st.secrets.get("models", {}).get("regular", "gemini-1.5-flash-latest")
    return genai.GenerativeModel(model_name)
@st.cache_resource
def get_pro_model():
    model_name = st.secrets.get("models", {}).get("pro", "gemini-pro")
    return genai.GenerativeModel(model_name)

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
    try: hoat_chat_goc = response_text.split("Output:")[1].strip()
    except IndexError: hoat_chat_goc = response_text.strip()
    if hoat_chat_goc == "INVALID" or not hoat_chat_goc: return f"❌ Lỗi: '{drug_name}' không được nhận dạng."
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

# --- 4. HÀM LOGIC TRUNG TÂM ---
def run_lookup(drug_name):
    try:
        is_pro = st.session_state.get("pro_access", False)
        final_result = get_drug_info(drug_name, is_pro_user=is_pro)
        if not final_result.startswith("❌ Lỗi:"):
            st.markdown(final_result)
            if drug_name not in st.session_state.history:
                st.session_state.history.insert(0, drug_name)
                if len(st.session_state.history) > 10: st.session_state.history.pop()
        else: st.error(final_result)
    except Exception as e:
        st.error("💥 Lỗi không xác định.")
        st.exception(e)

# --- 5. GIAO DIỆN VÀ LOGIC CHÍNH ---
st.set_page_config(page_title="Dược Điển AI", page_icon="💊")

# --- MÃ MỚI: TẢI CẤU HÌNH VÀ KHỞI TẠO TRÌNH XÁC THỰC ---
try:
    with open('config.yaml') as file:
        config = yaml.load(file, Loader=SafeLoader)
except FileNotFoundError:
    st.error("LỖI: Không tìm thấy file 'config.yaml'. Vui lòng tạo file này.")
    st.stop()

authenticator = stauth.Authenticate(
    config['credentials'],
    config['cookie']['name'],
    config['cookie']['key'],
    config['cookie']['expiry_days']
    # THAM SỐ 'pre_authorized' ĐÃ BỊ XÓA
)

# --- MÃ MỚI: HIỂN THỊ WIDGET ĐĂNG NHẬP / ĐĂNG KÝ ---
with st.sidebar:
    name, authentication_status, username = authenticator.login()
    if st.session_state["authentication_status"]:
        st.write(f'Chào mừng *{st.session_state["name"]}*')
        authenticator.logout('Đăng xuất')
    elif st.session_state["authentication_status"] is False:
        st.error('Tên đăng nhập/mật khẩu không chính xác')
    elif st.session_state["authentication_status"] is None:
        st.warning('Vui lòng nhập tên đăng nhập và mật khẩu')
        
# --- BẮT ĐẦU CẤU TRÚC PHÂN LUỒNG ---
st.title("Dược Điển AI 💊")
st.caption("Dự án được phát triển bởi group CÂCK và AI")


# --- LUỒNG 1: DÀNH CHO NGƯỜI DÙNG ĐÃ ĐĂNG NHẬP ---
if st.session_state["authentication_status"]:
    st.sidebar.success("Bạn đã đăng nhập.")
    st.sidebar.markdown("---")
    
    # --- Phần giao diện chính cho người dùng đã đăng nhập ---
    st.sidebar.header("Lịch sử tra cứu")
    if not st.session_state.history:
        st.sidebar.info("Chưa có thuốc nào được tra cứu.")
    else:
        for drug in st.session_state.history:
            if st.sidebar.button(drug, key=f"history_{drug}", use_container_width=True):
                run_lookup(drug)
    
    st.sidebar.markdown("---")
    
    with st.sidebar.container(border=True):
        st.write("**Bạn có ý tưởng để cải thiện ứng dụng?**")
        st.link_button( "Gửi phản hồi ngay!", url="https://forms.gle/M44GDS4hJ7LpY7b98", help="Mở form góp ý trong một tab mới" )
    
    st.sidebar.markdown("---")

    # Hiển thị tra cứu pro nếu có quyền
    st.sidebar.header("Truy cập Pro")
    if st.session_state.get("pro_access"):
        st.sidebar.success("Bạn đã có quyền truy cập Pro.")
    else:
        pro_code_input = st.sidebar.text_input("Nhập mã truy cập Pro:", type="password", help="Nhập mã của bạn và bấm nút Xác thực.")
        if st.sidebar.button("Xác thực"):
            is_valid, message = verify_code(pro_code_input)
            if is_valid:
                st.sidebar.success(message)
                st.rerun()
            else:
                st.sidebar.error(message)

    drug_name_input = st.text_input("Nhập tên thuốc (biệt dược hoặc hoạt chất):", key="main_input")
    lookup_button = st.button("Tra cứu")

    if lookup_button:
        if not drug_name_input:
            st.warning("Vui lòng nhập tên thuốc trước khi tra cứu.")
        else:
            run_lookup(drug_name_input)


# --- LUỒNG 2: DÀNH CHO KHÁCH (CHƯA ĐĂNG NHẬP) ---
else:
    st.info("Chào mừng bạn đến với Dược Điển AI. Vui lòng đăng nhập từ thanh bên để sử dụng đầy đủ tính năng.")
    
    # Chúng ta vẫn cho phép khách tra cứu cơ bản
    st.markdown("### Tra cứu cơ bản")
    drug_name_input_guest = st.text_input("Nhập tên thuốc (biệt dược hoặc hoạt chất):", key="guest_input")
    lookup_button_guest = st.button("Tra cứu", key="guest_button")

    if lookup_button_guest:
        if not drug_name_input_guest:
            st.warning("Vui lòng nhập tên thuốc trước khi tra cứu.")
        else:
            # Khách chỉ được dùng tính năng tra cứu không-pro
            run_lookup(drug_name_input_guest)
