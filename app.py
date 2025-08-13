import streamlit as st
import google.generativeai as genai

# --- 1. CẤU HÌNH VÀ PROMPT ---

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
10. **Liều dùng:** (Trình bày liều dùng cụ thể cho các chỉ định chính và các đối tượng đặc biệt nếu có, ví dụ: người lớn, trẻ em, người suy gan, người suy thận. Trình bày dưới dạng bảng hoặc gạch đầu dòng nếu có thể để dễ so sánh.)**

**QUY TẮC BẮT BUỘC:**
- Tuyệt đối KHÔNG được bịa đặt hay suy diễn thông tin.
- Nếu không tìm thấy dữ liệu cho mục nào, hãy ghi rõ: "Không có đủ dữ liệu đáng tin cậy."
- Luôn ưu tiên thông tin được chấp thuận bởi FDA.
"""

# --- 2. GIAO DIỆN NGƯỜI DÙNG ---

st.title("Dược Điển AI")
st.caption("Dự án được phát triển bởi group CÂCK và AI") # <--- THÊM DÒNG NÀY

drug_name = st.text_input("Nhập tên thuốc:")
lookup_button = st.button("Tra cứu")

# --- 3. LOGIC CỐT LÕI ---

if lookup_button:
    if not drug_name:
        st.warning("Vui lòng nhập tên thuốc trước khi tra cứu.")
    else:
        # Khối try-except đơn giản để bắt mọi lỗi
        try:
            with st.spinner("Dược sĩ AI đang tổng hợp thông tin..."):
                model = genai.GenerativeModel('gemini-2.5-flash')
                full_prompt = f"{PROMPT_GOC_RUT_GON}\n\nHãy tra cứu và trình bày thông tin cho thuốc sau đây: **{drug_name}**"
                response = model.generate_content(full_prompt)
                st.markdown(response.text)

        except Exception as e:
            # Hiển thị một thông báo lỗi chung chung cho MỌI VẤN ĐỀ
            st.error("Rất tiếc, đã có lỗi xảy ra trong quá trình tra cứu. Vui lòng thử lại sau ít phút.")
            # Dòng sau giúp chúng ta xem lỗi chi tiết là gì, nhưng người dùng không cần thấy
            st.exception(e)




