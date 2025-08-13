import streamlit as st
import google.generativeai as genai
import time
from google.api_core import exceptions

# ==============================================================================
# PHẦN CẤU HÌNH AI - SỬ DỤNG API KEY TỪ FILE SECRETS
# ==============================================================================
try:
    # Lấy API key từ file .streamlit/secrets.toml một cách an toàn
    genai.configure(api_key=st.secrets["GOOGLE_API_KEY"])
except FileNotFoundError:
    st.error("LỖI: Không tìm thấy file 'secrets.toml'. Vui lòng tạo file .streamlit/secrets.toml và thêm GOOGLE_API_KEY của bạn vào đó.")
    st.stop()
except KeyError:
    st.error("LỖI: Không tìm thấy 'GOOGLE_API_KEY' trong file 'secrets.toml'. Vui lòng kiểm tra lại file cấu hình.")
    st.stop()

# ==============================================================================
# PROMPT GỐC - BỘ NÃO CỦA DƯỢC ĐIỂN AI
# Đây là prompt đã được định nghĩa trong tài liệu dự án [cite: 59]
# ==============================================================================
PROMPT_GOC = """
Bạn là một Dược sĩ lâm sàng AI chuyên nghiệp và là chuyên gia trong việc tổng hợp thông tin y khoa. [cite: 61]
Nhiệm vụ của bạn là tra cứu và phân tích thông tin về một loại thuốc mà tôi cung cấp. [cite: 62]
Hãy sử dụng toàn bộ kiến thức đã được huấn luyện của bạn từ các nguồn dữ liệu y khoa uy tín trên thế giới như sách giáo khoa (Goodman & Gilman's, Katzung's), các cơ sở dữ liệu mở (openFDA, WHO), và các tạp chí khoa học hàng đầu (PubMed, The Lancet, NEJM). [cite: 63]

Khi tôi đưa tên một loại thuốc (có thể là tên gốc hoặc biệt dược), bạn PHẢI trình bày kết quả theo đúng cấu trúc 11 mục sau đây, sử dụng ngôn ngữ chuyên môn, chính xác và rõ ràng: [cite: 64]

1.  **Tên thuốc:** (Tên gốc và các tên biệt dược phổ biến) [cite: 65]
2.  **Nhóm thuốc:** [cite: 66]
3.  **Cơ chế:** [cite: 67]
4.  **Dược động học (ADME):** (Trình bày đủ từng mục A/D/M/E) [cite: 68]
5.  **Chỉ định:** [cite: 69]
6.  **Chống chỉ định:** [cite: 70]
7.  **Tương tác thuốc:** [cite: 71]
8.  **Tác dụng phụ:** [cite: 72]
9.  **Lưu ý lâm sàng & Theo dõi:** [cite: 73]
10. **Liều dùng:** [cite: 74]
11. **Nghiên cứu mới nhất:** (Tóm tắt 5-6 nghiên cứu nổi bật từ PubMed hoặc các tạp chí uy tín (Lancet, NEJM, Nature Reviews Drug Discovery, Pharmacological Reviews, Frontiers in Pharmacology, Drug Resistance Updates, Pharmacological Research) liên quan đến thuốc trong 1-2 năm gần nhất. Trình bày logic khoa học dễ hiểu) [cite: 75]

**QUY TẮC BẮT BUỘC:**
- Tuyệt đối KHÔNG được bịa đặt hay suy diễn thông tin. [cite: 77]
- Nếu không tìm thấy dữ liệu cho mục nào, hãy ghi rõ: "Không có đủ dữ liệu đáng tin cậy." [cite: 78]
- Luôn ưu tiên thông tin được chấp thuận bởi FDA. [cite: 79]
"""

# ==============================================================================
# (Phần mã xử lý lỗi và giao diện giữ nguyên như cũ)
# ==============================================================================
if 'button_disabled' not in st.session_state:
    st.session_state.button_disabled = False
if 'error_time' not in st.session_state:
    st.session_state.error_time = 0.0

st.title("Dược Điển AI - Phiên bản Thử nghiệm")
drug_name = st.text_input("Nhập tên thuốc (ví dụ: Atorvastatin, Paracetamol):")

COOLDOWN_SECONDS = 60
time_since_error = time.time() - st.session_state.error_time

if st.session_state.button_disabled and time_since_error < COOLDOWN_SECONDS:
    remaining_time = int(COOLDOWN_SECONDS - time_since_error)
    st.warning(f"💡 Lượng truy cập đang tạm thời quá tải. Vui lòng thử lại sau {remaining_time} giây.")
else:
    st.session_state.button_disabled = False

lookup_button = st.button("Tra cứu Thông Tin Thuốc", disabled=st.session_state.button_disabled)

# ==============================================================================
# PHẦN MÃ ĐƯỢC NÂNG CẤP - GỌI AI THẬT
# ==============================================================================
if lookup_button and drug_name:
    try:
        with st.spinner("Dược sĩ AI đang tổng hợp thông tin, vui lòng chờ... Điều này có thể mất một lúc."):
            # Thiết lập mô hình AI
            model = genai.GenerativeModel('gemini-1.5-pro')
            
            # Tạo câu lệnh hoàn chỉnh để gửi cho AI
            full_prompt = f"{PROMPT_GOC}\n\nHãy tra cứu và trình bày thông tin cho thuốc sau đây: **{drug_name}**"
            
            # Gửi yêu cầu đến AI và nhận kết quả
            response = model.generate_content(full_prompt)
            
            # Hiển thị kết quả ra màn hình
            st.markdown(response.text)

    except exceptions.ResourceExhausted as e:
        st.session_state.button_disabled = True
        st.session_state.error_time = time.time()
        st.error("Rất tiếc, đã có lỗi xảy ra do quá tải. Hệ thống sẽ tự động thử lại sau ít phút.")
        st.exception(e)
        st.rerun()
    except Exception as e:
        st.error("Đã có lỗi không xác định xảy ra. Vui lòng kiểm tra lại API Key và kết nối mạng.")
        st.exception(e)
