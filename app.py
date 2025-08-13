import streamlit as st
import google.generativeai as genai
import time
from google.api_core import exceptions

# --- 1. CẤU HÌNH VÀ HẰNG SỐ ---

# Cấu hình AI với API Key từ file secrets
try:
    genai.configure(api_key=st.secrets["GOOGLE_API_KEY"])
except (FileNotFoundError, KeyError):
    st.error("LỖI: Vui lòng tạo file .streamlit/secrets.toml và thêm `GOOGLE_API_KEY = 'KEY_CUA_BAN'` vào đó.")
    st.stop()

# Prompt gốc - Phiên bản rút gọn
PROMPT_GOC_RUT_GON = """
Bạn là một Dược sĩ lâm sàng AI chuyên nghiệp và là chuyên gia trong việc tổng hợp thông tin y khoa.
Nhiệm vụ của bạn là tra cứu và phân tích thông tin về một loại thuốc mà tôi cung cấp.
Hãy sử dụng toàn bộ kiến thức đã được huấn luyện của bạn từ các nguồn dữ liệu y khoa uy tín trên thế giới như sách giáo khoa (Goodman & Gilman's, Katzung's), các cơ sở dữ liệu mở (openFDA, WHO), và các tạp chí khoa học hàng đầu (PubMed, The Lancet, NEJM).

Khi tôi đưa tên một loại thuốc (có thể là tên gốc hoặc biệt dược), bạn PHẢI trình bày kết quả theo đúng cấu trúc 10 mục sau đây, sử dụng ngôn ngữ chuyên môn, chính xác và rõ ràng:

1.  **Tên thuốc:** (Tên gốc và các tên biệt dược phổ biến)
2.  **Nhóm thuốc:**
3.  **Cơ chế:**
4.  **Dược động học (ADME):** (Trình bày đủ từng mục A/D/M/E)
5.  **Chỉ định:**
6.  **Chống chỉ định:**
7.  **Tương tác thuốc:**
8.  **Tác dụng phụ:**
9.  **Lưu ý lâm sàng & Theo dõi:**
10. **Liều dùng:**

**QUY TẮC BẮT BUỘC:**
- Tuyệt đối KHÔNG được bịa đặt hay suy diễn thông tin.
- Nếu không tìm thấy dữ liệu cho mục nào, hãy ghi rõ: "Không có đủ dữ liệu đáng tin cậy."
- Luôn ưu tiên thông tin được chấp thuận bởi FDA.
"""

# Thời gian chờ khi gặp lỗi quá tải
COOLDOWN_SECONDS = 61 # Tăng thêm 1 giây để đảm bảo máy chủ sẵn sàng

# --- 2. QUẢN LÝ TRẠNG THÁI (SESSION STATE) ---

if 'button_disabled' not in st.session_state:
    st.session_state.button_disabled = False
if 'last_error_time' not in st.session_state:
    st.session_state.last_error_time = 0

# --- 3. GIAO DIỆN NGƯỜI DÙNG ---

st.title("Dược Điển AI - Phiên bản Thử nghiệm")
drug_name = st.text_input("Nhập tên thuốc (ví dụ: Atorvastatin, Paracetamol):")

# --- LOGIC KHÓA NÚT BẤM ---
# Kiểm tra xem có đang trong thời gian chờ không
elapsed_time = time.time() - st.session_state.last_error_time
if elapsed_time < COOLDOWN_SECONDS:
    st.session_state.button_disabled = True
    remaining_time = int(COOLDOWN_SECONDS - elapsed_time)
    st.warning(f"💡 Lượng truy cập đang tạm thời quá tải. Vui lòng thử lại sau {remaining_time} giây.")
else:
    st.session_state.button_disabled = False
    st.session_state.last_error_time = 0

# Nút bấm được điều khiển bởi session_state
lookup_button = st.button("Tra cứu Thông Tin Thuốc", disabled=st.session_state.button_disabled)


# --- 4. LOGIC XỬ LÝ CHÍNH ---

if lookup_button:
    if not drug_name:
        st.warning("Vui lòng nhập tên thuốc trước khi tra cứu.")
    else:
        try:
            with st.spinner("Dược sĩ AI đang tổng hợp thông tin, vui lòng chờ..."):
                model = genai.GenerativeModel('gemini-1.5-pro')
                full_prompt = f"{PROMPT_GOC_RUT_GON}\n\nHãy tra cứu và trình bày thông tin cho thuốc sau đây: **{drug_name}**"
                response = model.generate_content(full_prompt)
                st.markdown(response.text)

        except exceptions.ResourceExhausted:
            # Ghi lại thời điểm lỗi và chạy lại giao diện để khóa nút
            st.session_state.last_error_time = time.time()
            st.rerun()

        except Exception as e:
            st.error("Đã có lỗi không xác định xảy ra. Vui lòng kiểm tra lại.")
            st.exception(e)
