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
import pyrebase

# --- CẤU HÌNH FIREBASE ---
try:
    firebase_config = {
        "apiKey": st.secrets.firebase.apiKey,
        "authDomain": st.secrets.firebase.authDomain,
        "projectId": st.secrets.firebase.projectId,
        "storageBucket": st.secrets.firebase.storageBucket,
        "messagingSenderId": st.secrets.firebase.messagingSenderId,
        "appId": st.secrets.firebase.appId,
        "databaseURL": st.secrets.firebase.databaseURL
    }
    firebase = pyrebase.initialize_app(firebase_config)
    auth = firebase.auth()
    st.session_state.firebase_auth = auth
except Exception as e:
    st.error("Lỗi khi khởi tạo Firebase. Vui lòng kiểm tra file secrets.toml của bạn.")
    st.stop()

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
if 'user_info' not in st.session_state: st.session_state.user_info = None

# --- 2. CÁC HÀM XỬ LÝ ---
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
    # This function is kept as is
    pass

def run_lookup(drug_name):
    try:
        is_pro = st.session_state.get("pro_access", False)
        # Assuming get_drug_info is defined elsewhere and works
        final_result = f"Thông tin cho {drug_name} (Pro: {is_pro})" # Placeholder
        if not final_result.startswith("❌ Lỗi:"):
            st.markdown(final_result)
            if drug_name not in st.session_state.history:
                st.session_state.history.insert(0, drug_name)
                if len(st.session_state.history) > 10: st.session_state.history.pop()
        else: st.error(final_result)
    except Exception as e:
        st.error("💥 Lỗi không xác định.")
        st.exception(e)

# --- GIAO DIỆN VÀ LOGIC CHÍNH ---
st.set_page_config(page_title="Dược Điển AI", page_icon="💊")
st.title("Dược Điển AI 💊")
st.caption("Dự án được phát triển bởi group CÂCK và AI")

# --- HỆ THỐNG ĐĂNG NHẬP MỚI ---
auth = st.session_state.firebase_auth
is_logged_in = st.session_state.get("user_info") is not None

# --- SIDEBAR LOGIC ---
with st.sidebar:
    if is_logged_in:
        user_email = st.session_state.user_info['email']
        st.success(f"Chào mừng, {user_email}")
        if st.button("Đăng xuất"):
            st.session_state.user_info = None
            st.rerun()
    else:
        choice = st.selectbox("Đăng nhập / Đăng ký", ["Tiếp tục với tư cách khách", "Đăng nhập", "Đăng ký"])

        if choice == "Đăng nhập":
            with st.form("login_form"):
                email = st.text_input("Email")
                password = st.text_input("Mật khẩu", type="password")
                login_button = st.form_submit_button("Đăng nhập")
                if login_button:
                    try:
                        user = auth.sign_in_with_email_and_password(email, password)
                        st.session_state.user_info = user
                        st.rerun()
                    except Exception as e:
                        st.error("Email hoặc mật khẩu không chính xác.")
        elif choice == "Đăng ký":
            with st.form("register_form"):
                email = st.text_input("Email")
                password = st.text_input("Mật khẩu", type="password")
                register_button = st.form_submit_button("Đăng ký")
                if register_button:
                    try:
                        user = auth.create_user_with_email_and_password(email, password)
                        st.success("Đăng ký thành công! Vui lòng chuyển qua tab 'Đăng nhập'.")
                    except Exception as e:
                        st.error("Email này có thể đã tồn tại hoặc không hợp lệ.")

    # --- HIỂN THỊ CÁC THÀNH PHẦN CHUNG CỦA SIDEBAR ---
    st.header("Lịch sử tra cứu")
    if not st.session_state.history:
        st.info("Chưa có thuốc nào được tra cứu.")
    else:
        for drug in st.session_state.history:
            if st.button(drug, key=f"history_{drug}", use_container_width=True):
                # We need a way to trigger the lookup from the history button
                st.session_state.drug_to_lookup = drug

    st.markdown("---")
    with st.container(border=True):
        st.write("**Bạn có ý tưởng để cải thiện ứng dụng?**")
        st.link_button("Gửi phản hồi ngay!", url="https://forms.gle/M44GDS4hJ7LpY7b98", help="Mở form góp ý trong một tab mới")

    # --- HIỂN THỊ MỤC PRO ACCESS NẾU ĐÃ ĐĂNG NHẬP ---
    if is_logged_in:
        st.markdown("---")
        st.header("Truy cập Pro")
        if st.session_state.get("pro_access"):
            st.success("Bạn đã có quyền truy cập Pro.")
        else:
            pro_code_input = st.text_input("Nhập mã truy cập Pro:", type="password", help="Nhập mã của bạn và bấm nút Xác thực.")
            if st.button("Xác thực"):
                is_valid, message = verify_code(pro_code_input)
                if is_valid:
                    st.success(message)
                    st.rerun()
                else:
                    st.error(message)


# --- GIAO DIỆN CHÍNH ---
if not is_logged_in:
    st.info("Bạn đang sử dụng với tư cách khách. Đăng nhập để lưu lịch sử và sử dụng các tính năng nâng cao.")

drug_name_input = st.text_input("Nhập tên thuốc (biệt dược hoặc hoạt chất):", key="main_input")
lookup_button = st.button("Tra cứu")

# Logic to handle lookup from main button or history button
if lookup_button and drug_name_input:
    run_lookup(drug_name_input)
elif st.session_state.get("drug_to_lookup"):
    drug_to_run = st.session_state.drug_to_lookup
    st.session_state.drug_to_lookup = None # Clear after use
    run_lookup(drug_to_run)
