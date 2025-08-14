import streamlit as st
import google.generativeai as genai
import google.generativeai.types as genai_types
import google.api_core.exceptions as ga_ex
import datetime
import pandas as pd
import gspread
from gspread_dataframe import get_as_dataframe
from gspread.exceptions import SpreadsheetNotFound, WorksheetNotFound

# --- KIỂM TRA TRẠNG THÁI BẢO TRÌ ---
is_maintenance = st.secrets.get("maintenance_mode", False) 
if is_maintenance:
    st.set_page_config(page_title="Bảo trì", page_icon="🛠️")
    st.title("🛠️ Dược Điển AI đang được bảo trì")
    message = st.secrets.get("maintenance_message", "Ứng dụng đang được cập nhật. Vui lòng quay lại sau.")
    st.info(message)
    st.stop()

# --- 1. KHỞI TẠO TRẠNG THÁI PHIÊN (SESSION STATE) ---
if 'history' not in st.session_state:
    st.session_state.history = []
if 'pro_access' not in st.session_state:
    st.session_state.pro_access = False

# --- 2. CẤU HÌNH VÀ TẢI PROMPTS ---
def load_prompt(file_path):
    try:
        with open(file_path, "r", encoding="utf-8") as f:
            return f.read()
    except FileNotFoundError:
        st.error(f"LỖI: Không tìm thấy file prompt tại '{file_path}'. Vui lòng đảm bảo file này tồn tại trong cùng thư mục với app.py.")
        st.stop()

try:
    genai.configure(api_key=st.secrets["GOOGLE_API_KEY"])
except (FileNotFoundError, KeyError):
    st.error("LỖI: Vui lòng tạo file .streamlit/secrets.toml và thêm `GOOGLE_API_KEY = 'KEY_CUA_BAN'` vào đó.")
    st.stop()

PROMPT_NHAN_DIEN = load_prompt("prompt_nhandien.txt")
PROMPT_GOC_RUT_GON = load_prompt("prompt_goc_rutgon.txt")

# --- 3. CÁC HÀM XỬ LÝ ---

# --- HÀM XỬ LÝ MÃ TRUY CẬP (PHIÊN BẢN GSPREAD DEBUG) ---
@st.cache_data(ttl=600)
def get_access_codes_df():
    """Kết nối tới Google Sheets bằng gspread và lấy dữ liệu với các bước debug."""
    try:
        st.info("DEBUG: Bắt đầu kết nối tới Google Sheets...")
        credentials = st.secrets.connections.gsheets.credentials
        gspread_client = gspread.service_account_from_dict(credentials)
        st.info("DEBUG: Xác thực với Google thành công.")
        
        spreadsheet_name = st.secrets.connections.gsheets.spreadsheet
        st.info(f"DEBUG: Đang thử mở spreadsheet có tên: '{spreadsheet_name}'...")
        spreadsheet = gspread_client.open(spreadsheet_name)
        st.info("DEBUG: Mở spreadsheet thành công. Đang đọc worksheet...")
        
        worksheet = spreadsheet.sheet1
        st.info(f"DEBUG: Đang đọc từ worksheet có tên: '{worksheet.title}'...")
        codes_df = get_as_dataframe(worksheet)
        st.info("DEBUG: Đọc dữ liệu thành công!")
        return codes_df

    except SpreadsheetNotFound:
        st.error(f"LỖI DEBUG: Không tìm thấy Spreadsheet có tên '{st.secrets.connections.gsheets.spreadsheet}'. Vui lòng kiểm tra lại tên trong secrets và đảm bảo bạn đã chia sẻ file cho email: {st.secrets.connections.gsheets.credentials.client_email}")
        return pd.DataFrame()
    except WorksheetNotFound:
        st.error("LỖI DEBUG: Không tìm thấy worksheet có tên 'Sheet1'. Có thể bạn đã đổi tên trang tính đầu tiên. Vui lòng đổi tên nó lại thành 'Sheet1'.")
        return pd.DataFrame()
    except Exception as e:
        st.error("LỖI DEBUG không xác định. Thông tin lỗi chi tiết:")
        st.exception(e) # In ra toàn bộ lỗi để chúng ta xem
        return pd.DataFrame()


def verify_code(user_code):
    """Kiểm tra mã người dùng nhập với dữ liệu trên Google Sheet."""
    if not user_code:
        return False, "Vui lòng nhập mã truy cập."
    
    codes_df = get_access_codes_df()
    if codes_df.empty:
        return False, "Không thể tải dữ liệu mã truy cập. Vui lòng kiểm tra các thông báo lỗi DEBUG ở trên."

    codes_df.dropna(subset=['code'], inplace=True)
    codes_df['code'] = codes_df['code'].astype(str)
    
    matched_code_series = codes_df[codes_df['code'] == user_code]

    if matched_code_series.empty:
        return False, "Mã không hợp lệ hoặc không tìm thấy."

    code_info = matched_code_series.iloc[0]
    code_type = code_info['type']
    
    if code_type == 'permanent':
        st.session_state.pro_access = True
        return True, f"Xác thực thành công! Chào mừng {code_info.get('owner', 'Pro User')}."

    if code_type == 'temporary':
        try:
            created_date = pd.to_datetime(code_info['created_at']).date()
            today = datetime.date.today()
            days_passed = (today - created_date).days
            
            if 0 <= days_passed <= 7:
                st.session_state.pro_access = True
                return True, f"Xác thực thành công! Mã của bạn còn hiệu lực {7 - days_passed} ngày."
            else:
                return False, "Mã tạm thời của bạn đã hết hạn."
        except Exception:
            return False, "Lỗi định dạng ngày tháng trong file Google Sheet. Vui lòng kiểm tra lại cột 'created_at'."
    
    return False, "Loại mã không xác định."

# --- HÀM XỬ LÝ DƯỢC ĐIỂN ---
# ... (Các hàm còn lại không thay đổi) ...
@st.cache_resource
def get_model():
    return genai.GenerativeModel('gemini-2.5-flash-lite')

@st.cache_data(ttl="6h")
def get_drug_info(drug_name, is_pro_user=False):
    model = get_model()
    prompt_nhan_dien_final = PROMPT_NHAN_DIEN.format(drug_name=drug_name)
    response_nhan_dien = model.generate_content(prompt_nhan_dien_final)
    
    response_text = response_nhan_dien.text
    try:
        hoat_chat_goc = response_text.split("Output:")[1].strip()
    except IndexError:
        hoat_chat_goc = response_text.strip()

    if hoat_chat_goc == "INVALID" or not hoat_chat_goc:
        return f"❌ Lỗi: '{drug_name}' không được nhận dạng là một tên thuốc hợp lệ."

    generation_config = {
        "max_output_tokens": 8192,
        "temperature": 0.6,
    }
    full_prompt = f"{PROMPT_GOC_RUT_GON}\n\nHãy tra cứu và trình bày thông tin cho thuốc sau đây: **{hoat_chat_goc}**"
    response_phan_tich = model.generate_content(full_prompt, generation_config=generation_config)
    final_response = f"✅ Hoạt chất đã nhận diện: **{hoat_chat_goc}**\n\n---\n\n{response_phan_tich.text}"
    
    if is_pro_user:
        final_response += "\n\n---\n\n**11. Nghiên cứu lâm sàng gần đây:** (Tính năng Pro. Sắp ra mắt...)"
        
    return final_response

# --- 4. HÀM LOGIC TRUNG TÂM ---
def run_lookup(drug_name):
    try:
        with st.spinner(f"Đang tra cứu '{drug_name}'..."):
            is_pro = st.session_state.get("pro_access", False)
            final_result = get_drug_info(drug_name, is_pro_user=is_pro)
        if not final_result.startswith("❌ Lỗi:"):
            st.markdown(final_result)
            if drug_name not in st.session_state.history:
                st.session_state.history.insert(0, drug_name)
                if len(st.session_state.history) > 10:
                     st.session_state.history.pop()
        else:
            st.error(final_result)
    except Exception as e:
        st.error("💥 Lỗi không xác định: Một sự cố không mong muốn đã xảy ra.")
        st.exception(e)

# --- 5. GIAO DIỆN VÀ LOGIC CHÍNH ---
st.set_page_config(page_title="Dược Điển AI", page_icon="💊")
st.title("Dược Điển AI 💊")
st.caption("Dự án được phát triển bởi group CÂCK và AI")

# --- Sidebar ---
st.sidebar.header("Lịch sử tra cứu")
# ... (Phần này giữ nguyên)
for drug in st.session_state.history:
    if st.sidebar.button(drug, key=f"history_{drug}", use_container_width=True):
        run_lookup(drug)

st.sidebar.markdown("---")
with st.sidebar.container(border=True):
    st.write("**Bạn có ý tưởng để cải thiện ứng dụng?**")
    st.link_button( "Gửi phản hồi ngay!", url="https://forms.gle/M44GDS4hJ7LpY7b98", help="Mở form góp ý trong một tab mới" )
st.sidebar.markdown("---")

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

# --- Main page ---
drug_name_input = st.text_input("Nhập tên thuốc (biệt dược hoặc hoạt chất):", key="main_input")
lookup_button = st.button("Tra cứu")

if lookup_button:
    if not drug_name_input:
        st.warning("Vui lòng nhập tên thuốc trước khi tra cứu.")
    else:
        run_lookup(drug_name_input)
