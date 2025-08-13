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

# Prompt gốc - Phiên bản rút gọn để giảm tải
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
COOLDOWN_SECONDS = 60


# --- 2. GIAO DIỆN NGƯỜI DÙNG ---

st.title("Dược Điển AI - Phiên bản Thử nghiệm")
drug_name = st.text_input("Nhập tên thuốc (ví dụ: Atorvastatin, Paracetamol):")

# Chỉ có MỘT nút bấm duy nhất trong toàn bộ ứng dụng
lookup_button = st.button("Tra cứu Thông Tin Thuốc")


# --- 3. LOGIC XỬ LÝ CHÍNH ---

# Toàn bộ logic sẽ chỉ chạy khi người dùng bấm nút
if lookup_button:
    # Kiểm tra xem người dùng đã nhập tên thuốc chưa
    if not drug_name:
        st.warning("Vui lòng nhập tên thuốc trước khi tra cứu.")
    else:
        # Bắt đầu xử lý chính khi đã có tên thuốc
        try:
            with st.spinner("Dược sĩ AI đang tổng hợp thông tin, vui lòng chờ..."):
                # Thiết lập mô hình AI
                model = genai.GenerativeModel('gemini-1.5-pro')

                # Tạo câu lệnh hoàn chỉnh để gửi cho AI (sử dụng prompt rút gọn)
                full_prompt = f"{PROMPT_GOC_RUT_GON}\n\nHãy tra cứu và trình bày thông tin cho thuốc sau đây: **{drug_name}**"

                # Gửi yêu cầu đến AI và nhận kết quả
                response = model.generate_content(full_prompt)

                # Hiển thị kết quả ra màn hình
                st.markdown(response.text)

        except exceptions.ResourceExhausted:
            # Xử lý lỗi quá tải với đồng hồ đếm ngược
            placeholder = st.empty()
            for i in range(COOLDOWN_SECONDS, 0, -1):
                placeholder.warning(f"💡 Lượng truy cập đang tạm thời quá tải. Vui lòng thử lại sau {i} giây.")
                time.sleep(1) # Chờ 1 giây
            placeholder.empty() # Xóa thông báo khi đếm ngược xong

        except Exception as e:
            # Xử lý các lỗi không xác định khác
            st.error("Đã có lỗi không xác định xảy ra. Vui lòng kiểm tra lại.")
            st.exception(e) # In ra lỗi chi tiết để chúng ta gỡ rối
