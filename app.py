import streamlit as st
import google.generativeai as genai
import google.generativeai.types as genai_types
import google.api_core.exceptions as ga_ex

# --- KIỂM TRA TRẠNG THÁI BẢO TRÌ ---
# Đọc "công tắc" từ secrets. Mặc định là False (không bảo trì) nếu không tìm thấy.
is_maintenance = st.secrets.get("maintenance_mode", False) 

if is_maintenance:
    st.set_page_config(page_title="Bảo trì", page_icon="🛠️")
    st.title("🛠️ Dược Điển AI đang được bảo trì")
    # Lấy thông báo bảo trì từ secrets, nếu không có thì dùng thông báo mặc định.
    message = st.secrets.get("maintenance_message", "Ứng dụng đang được cập nhật. Vui lòng quay lại sau.")
    st.info(message)
    st.stop() # Dừng toàn bộ phần còn lại của ứng dụng không cho chạy.

# --- 1. KHỞI TẠO TRẠNG THÁI PHIÊN (SESSION STATE) ---
if 'history' not in st.session_state:
    st.session_state.history = []

# --- 2. CẤU HÌNH VÀ TẢI PROMPTS ---
def load_prompt(file_path):
    """Hàm này đọc nội dung từ một file và trả về dưới dạng chuỗi."""
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

# Tải các prompt từ file bên ngoài
PROMPT_NHAN_DIEN = load_prompt("prompt_nhandien.txt")
PROMPT_GOC_RUT_GON = load_prompt("prompt_goc_rutgon.txt")


# --- 3. CÁC HÀM XỬ LÝ (Cache) ---
@st.cache_resource
def get_model():
    return genai.GenerativeModel('gemini-2.5-flash-lite')

@st.cache_data(ttl="6h")
def get_drug_info(drug_name):
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
    response_phan_tich = model.generate_content(
        full_prompt,
        generation_config=generation_config
    )
    final_response = f"✅ Hoạt chất đã nhận diện: **{hoat_chat_goc}**\n\n---\n\n{response_phan_tich.text}"
    return final_response
    
# --- 4. HÀM LOGIC TRUNG TÂM ---
def run_lookup(drug_name):
    try:
        with st.spinner(f"Đang tra cứu '{drug_name}'..."):
            final_result = get_drug_info(drug_name)
        if not final_result.startswith("❌ Lỗi:"):
            st.markdown(final_result)
            if drug_name not in st.session_state.history:
                st.session_state.history.insert(0, drug_name)
                if len(st.session_state.history) > 10:
                     st.session_state.history.pop()
        else:
            st.error(final_result)
    except ga_ex.PermissionDenied as e:
        st.error("🚫 Lỗi Xác Thực: Google API Key của bạn không hợp lệ hoặc đã bị vô hiệu hóa.")
    except ga_ex.ResourceExhausted as e:
        st.error("🚦 Đã đạt giới hạn: Bạn đã gửi quá nhiều yêu cầu trong một thời gian ngắn.")
    except ValueError as e:
        if "safety setting" in str(e):
            st.error("🔒 Nội dung bị chặn: Yêu cầu của bạn có thể đã vi phạm chính sách an toàn.")
        else:
            st.error(f"Lỗi Dữ Liệu: Có vấn đề với dữ liệu đầu vào hoặc đầu ra.")
    except ga_ex.GoogleAPICallError as e:
        st.error("🌐 Lỗi Kết Nối: Máy chủ Google AI đang gặp sự cố tạm thời. Vui lòng thử lại sau ít phút.")
    except Exception as e:
        st.error("💥 Lỗi không xác định: Một sự cố không mong muốn đã xảy ra.")
        st.exception(e)

# --- 5. GIAO DIỆN VÀ LOGIC CHÍNH ---
st.set_page_config(page_title="Dược Điển AI", page_icon="💊")
st.title("Dược Điển AI 💊")
st.caption("Dự án được phát triển bởi group CÂCK và AI")

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
    st.link_button(
        "Gửi phản hồi ngay!",
        url="https://forms.gle/M44GDS4hJ7LpY7b98",
        help="Mở form góp ý trong một tab mới"
    )

drug_name_input = st.text_input("Nhập tên thuốc (biệt dược hoặc hoạt chất):", key="main_input")
lookup_button = st.button("Tra cứu")

if lookup_button:
    if not drug_name_input:
        st.warning("Vui lòng nhập tên thuốc trước khi tra cứu.")
    else:
        run_lookup(drug_name_input)
