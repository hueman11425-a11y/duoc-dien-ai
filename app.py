import streamlit as st
import google.generativeai as genai

# Import các file chúng ta đã tách
import auth
import utils

# --- KHỞI TẠO CÁC DỊCH VỤ ---
# Khởi tạo Firebase Authentication
firebase_auth = auth.initialize_firebase()
if not firebase_auth:
    st.stop() # Dừng ứng dụng nếu không kết nối được Firebase

# Cấu hình Google AI
try:
    genai.configure(api_key=st.secrets["GOOGLE_API_KEY"])
except (FileNotFoundError, KeyError):
    st.error("LỖI: Vui lòng cấu hình GOOGLE_API_KEY trong secrets.toml.")
    st.stop()

# --- KIỂM TRA TRẠNG THÁI BẢO TRÌ ---
is_maintenance = st.secrets.get("maintenance_mode", False)
if is_maintenance:
    st.set_page_config(page_title="Bảo trì", page_icon="🛠️")
    st.title("🛠️ Dược Điển AI đang được bảo trì")
    message = st.secrets.get("maintenance_message", "Ứng dụng đang được cập nhật. Vui lòng quay lại sau.")
    st.info(message)
    st.stop()

# --- KHỞI TẠO TRẠNG THÁI PHIÊN ---
if 'history' not in st.session_state: st.session_state.history = []
if 'pro_access' not in st.session_state: st.session_state.pro_access = False

# --- HÀM LOGIC TRUNG TÂM ---
def run_lookup(drug_name):
    try:
        is_pro = st.session_state.get("pro_access", False)
        final_result = utils.get_drug_info(drug_name, is_pro_user=is_pro)
        if not final_result.startswith("❌ Lỗi:"):
            st.markdown(final_result)
            # Chỉ lưu vào lịch sử tạm thời nếu chưa đăng nhập
            if st.session_state.get("user_info") is None:
                if drug_name not in st.session_state.history:
                    st.session_state.history.insert(0, drug_name)
                    if len(st.session_state.history) > 10:
                        st.session_state.history.pop()
        else:
            st.error(final_result)
    except Exception as e:
        st.error("💥 Lỗi không xác định.")
        st.exception(e)

# --- BẮT ĐẦU GIAO DIỆN ---
st.set_page_config(page_title="Dược Điển AI", page_icon="💊")

# Hiển thị form đăng nhập và lấy trạng thái
is_logged_in = auth.display_auth_forms(firebase_auth)

# --- GIAO DIỆN CHÍNH ---
st.title("Dược Điển AI 💊")
st.caption("Dự án được phát triển bởi group CÂCK và AI")

# --- SIDEBAR (phần còn lại) ---
with st.sidebar:
    # Lịch sử và phản hồi hiển thị cho tất cả mọi người
    st.header("Lịch sử tra cứu")
    if not st.session_state.history:
        st.info("Chưa có thuốc nào được tra cứu.")
    else:
        for drug in st.session_state.history:
            if st.button(drug, key=f"history_{drug}", use_container_width=True):
                run_lookup(drug)

    st.markdown("---")
    with st.container(border=True):
        st.write("**Bạn có ý tưởng để cải thiện ứng dụng?**")
        st.link_button("Gửi phản hồi ngay!", url="https://forms.gle/M44GDS4hJ7LpY7b98", help="Mở form góp ý trong một tab mới")

    # Mục Pro chỉ hiển thị khi đã đăng nhập
    if is_logged_in:
        st.markdown("---")
        st.header("Truy cập Pro")
        if st.session_state.get("pro_access"):
            st.success("Bạn đã có quyền truy cập Pro.")
        else:
            pro_code_input = st.text_input("Nhập mã truy cập Pro:", type="password", help="Nhập mã của bạn và bấm nút Xác thực.")
            if st.button("Xác thực"):
                is_valid, message = utils.verify_code(pro_code_input)
                if is_valid:
                    st.success(message)
                    st.rerun()
                else:
                    st.error(message)

# --- KHUNG NHẬP LIỆU CHÍNH ---
if not is_logged_in:
    st.info("Bạn đang sử dụng với tư cách khách. Đăng nhập để lưu lịch sử và sử dụng các tính năng nâng cao.")

drug_name_input = st.text_input("Nhập tên thuốc (biệt dược hoặc hoạt chất):", key="main_input")
lookup_button = st.button("Tra cứu")

if lookup_button:
    if not drug_name_input:
        st.warning("Vui lòng nhập tên thuốc trước khi tra cứu.")
    else:
        run_lookup(drug_name_input)
