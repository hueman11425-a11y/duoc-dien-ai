import streamlit as st
import google.generativeai as genai
# --- IMPORT MỚI ĐỂ XỬ LÝ LỖI ---
import google.generativeai.types as genai_types
import google.api_core.exceptions as ga_ex

# --- 1. CẤU HÌNH VÀ PROMPTS ---
# (Không thay đổi)
try:
    genai.configure(api_key=st.secrets["GOOGLE_API_KEY"])
except (FileNotFoundError, KeyError):
    st.error("LỖI: Vui lòng tạo file .streamlit/secrets.toml và thêm `GOOGLE_API_KEY = 'KEY_CUA_BAN'` vào đó.")
    st.stop()
    
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


# --- 2. CÁC HÀM XỬ LÝ (Cache) ---
# (Không thay đổi)
@st.cache_resource
def get_model():
    print("--- Khởi tạo model AI ---")
    return genai.GenerativeModel('gemini-2.5-flash-lite')

@st.cache_data
def get_drug_info(drug_name):
    print(f"--- Thực hiện tra cứu API cho: {drug_name} ---")
    model = get_model()

    prompt_nhan_dien_final = PROMPT_NHAN_DIEN.format(drug_name=drug_name)
    response_nhan_dien = model.generate_content(prompt_nhan_dien_final)
    hoat_chat_goc = response_nhan_dien.text.strip()

    if not hoat_chat_goc or "không tìm thấy" in hoat_chat_goc.lower():
        return f"Lỗi: Không thể xác định được hoạt chất cho '{drug_name}'."

    full_prompt = f"{PROMPT_GOC_RUT_GON}\n\nHãy tra cứu và trình bày thông tin cho thuốc sau đây: **{hoat_chat_goc}**"
    response_phan_tich = model.generate_content(full_prompt)
    
    final_response = f"✅ Đã nhận diện hoạt chất: **{hoat_chat_goc}**\n\n---\n\n{response_phan_tich.text}"
    return final_response


# --- 3. GIAO DIỆN VÀ LOGIC CHÍNH ---

st.title("Dược Điển AI")
st.caption("Dự án được phát triển bởi group CÂCK và AI")

drug_name_input = st.text_input("Nhập tên thuốc (biệt dược hoặc hoạt chất):")
lookup_button = st.button("Tra cứu")

if lookup_button:
    if not drug_name_input:
        st.warning("Vui lòng nhập tên thuốc trước khi tra cứu.")
    else:
        # --- KHỐI TRY-EXCEPT ĐÃ ĐƯỢC NÂNG CẤP ---
        try:
            with st.spinner("Dược sĩ AI đang làm việc, vui lòng chờ..."):
                final_result = get_drug_info(drug_name_input)
            st.markdown(final_result)
        
        # Bắt lỗi cụ thể
        except ga_ex.PermissionDenied as e:
            st.error("🚫 Lỗi Xác Thực: Google API Key của bạn không hợp lệ hoặc đã bị vô hiệu hóa. Vui lòng kiểm tra lại trong file `.streamlit/secrets.toml`.")
            st.exception(e)

        except ga_ex.ResourceExhausted as e:
            st.error("🚦 Đã đạt giới hạn: Bạn đã gửi quá nhiều yêu cầu trong một thời gian ngắn. Vui lòng chờ vài phút rồi thử lại.")
            st.exception(e)
        
        except ValueError as e:
            # Lỗi này thường xảy ra khi nội dung bị chặn do chính sách an toàn
            if "safety setting" in str(e):
                st.error("🔒 Nội dung bị chặn: Yêu cầu của bạn có thể đã vi phạm chính sách an toàn của Google. Vui lòng thử lại với một tên thuốc khác.")
                st.exception(e)
            else:
                st.error(f"Lỗi Dữ Liệu: Có vấn đề với dữ liệu đầu vào hoặc đầu ra. Chi tiết: {e}")
                st.exception(e)

        except ga_ex.GoogleAPICallError as e:
            st.error("🌐 Lỗi Kết Nối: Không thể kết nối đến máy chủ của Google AI. Vui lòng kiểm tra lại kết nối mạng của bạn.")
            st.exception(e)

        # Bắt tất cả các lỗi còn lại
        except Exception as e:
            st.error("💥 Lỗi không xác định: Một sự cố không mong muốn đã xảy ra.")
            st.exception(e)

