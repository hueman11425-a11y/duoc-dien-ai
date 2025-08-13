import streamlit as st
import google.generativeai as genai

# --- 1. CẤU HÌNH VÀ PROMPTS ---

# Cấu hình AI với API Key từ file secrets
try:
    genai.configure(api_key=st.secrets["GOOGLE_API_KEY"])
except (FileNotFoundError, KeyError):
    st.error("LỖI: Vui lòng tạo file .streamlit/secrets.toml và thêm `GOOGLE_API_KEY = 'KEY_CUA_BAN'` vào đó.")
    st.stop()

# Prompt MỚI: Chỉ dùng để nhận diện hoạt chất từ biệt dược
PROMPT_NHAN_DIEN = """Từ tên thuốc sau đây, hãy trích xuất (các) hoạt chất gốc. 
Chỉ trả về tên (các) hoạt chất, phân cách bằng dấu phẩy nếu có nhiều hoạt chất. 
TUYỆT ĐỐI không giải thích, không thêm bất kỳ từ nào khác.

Ví dụ 1:
Input: Amoxicillin 500mg
Output: Amoxicillin

Ví dụ 2:
Input: Troysar AM
Output: Losartan, Amlodipine

Input: {drug_name}
Output:"""


# Prompt gốc - Dùng để phân tích chi tiết sau khi đã có tên hoạt chất
PROMPT_GOC_RUT_GON = """
Bạn là một Dược sĩ lâm sàng AI chuyên nghiệp và là chuyên gia trong việc tổng hợp thông tin y khoa.
Nhiệm vụ của bạn là tra cứu và phân tích thông tin về một loại thuốc mà tôi cung cấp.
Hãy sử dụng toàn bộ kiến thức đã được huấn luyện của bạn từ các nguồn dữ liệu y khoa uy tín trên thế giới như sách giáo khoa (Goodman & Gilman's, Katzung's), các cơ sở dữ liệu mở (openFDA, WHO), và các tạp chí khoa học hàng đầu (PubMed, The Lancet, NEJM).

Khi tôi đưa tên một loại thuốc (luôn là tên gốc/hoạt chất), bạn PHẢI trình bày kết quả theo đúng cấu trúc 10 mục sau đây, sử dụng ngôn ngữ chuyên môn, chính xác và rõ ràng:

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

st.title("Dược Điển AI (Bản nâng cấp)")
st.caption("Dự án được phát triển bởi group CÂCK và AI")

drug_name_input = st.text_input("Nhập tên thuốc (biệt dược hoặc hoạt chất):")
lookup_button = st.button("Tra cứu")

# --- 3. LOGIC CỐT LÕI (ĐÃ NÂNG CẤP) ---

if lookup_button:
    if not drug_name_input:
        st.warning("Vui lòng nhập tên thuốc trước khi tra cứu.")
    else:
        try:
            model = genai.GenerativeModel('gemini-2.5-flash-lite')

            # --- BƯỚC 1: NHẬN DIỆN HOẠT CHẤT ---
            with st.spinner(f"Đang nhận diện hoạt chất trong '{drug_name_input}'..."):
                prompt_nhan_dien_final = PROMPT_NHAN_DIEN.format(drug_name=drug_name_input)
                response_nhan_dien = model.generate_content(prompt_nhan_dien_final)
                # Dọn dẹp output, chỉ lấy text và xóa khoảng trắng thừa
                hoat_chat_goc = response_nhan_dien.text.strip()

            if not hoat_chat_goc or "không tìm thấy" in hoat_chat_goc.lower():
                 st.error(f"Không thể xác định được hoạt chất cho '{drug_name_input}'. Vui lòng thử lại với tên khác.")
            else:
                st.info(f"✅ Đã nhận diện hoạt chất: **{hoat_chat_goc}**")

                # --- BƯỚC 2: PHÂN TÍCH CHI TIẾT HOẠT CHẤT ĐÃ TÌM ĐƯỢC ---
                with st.spinner(f"Dược sĩ AI đang tổng hợp thông tin chi tiết về '{hoat_chat_goc}'..."):
                    full_prompt = f"{PROMPT_GOC_RUT_GON}\n\nHãy tra cứu và trình bày thông tin cho thuốc sau đây: **{hoat_chat_goc}**"
                    response_phan_tich = model.generate_content(full_prompt)
                    st.markdown(response_phan_tich.text)

        except Exception as e:
            st.error("Rất tiếc, đã có lỗi xảy ra trong quá trình tra cứu. Vui lòng thử lại sau ít phút.")
            st.exception(e) # Dòng này giúp bạn thấy lỗi chi tiết khi lập trình

