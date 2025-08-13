import streamlit as st
import google.generativeai as genai
import google.generativeai.types as genai_types
import google.api_core.exceptions as ga_ex

# --- 1. CẤU HÌNH VÀ PROMPTS ---
try:
    genai.configure(api_key=st.secrets["GOOGLE_API_KEY"])
except (FileNotFoundError, KeyError):
    st.error("LỖI: Vui lòng tạo file .streamlit/secrets.toml và thêm `GOOGLE_API_KEY = 'KEY_CUA_BAN'` vào đó.")
    st.stop()
    
PROMPT_NHAN_DIEN = """Từ tên thuốc sau đây, hãy trích xuất (các) hoạt chất gốc. 
Chỉ trả về tên (các) hoạt chất, phân cách bằng dấu phẩy nếu có nhiều hoạt chất. 
Nếu Input không giống một tên thuốc, hãy trả về duy nhất từ "LỖI".
TUYỆT ĐỐI không giải thích, không thêm bất kỳ từ nào khác.
Input: {drug_name}
Output:"""

# --- PROMPT MỚI ĐỂ XÁC THỰC ---
PROMPT_XAC_THUC = """Bạn là chuyên gia xác thực tên thuốc. 
Hãy xem xét liệu '{original_input}' có phải là một tên biệt dược, tên gốc, hoặc một cách gọi phổ biến của thuốc chứa hoạt chất '{identified_ingredient}' hay không.
Chỉ trả lời bằng một từ duy nhất: CÓ hoặc KHÔNG.

Ví dụ 1:
original_input: Viagra
identified_ingredient: Sildenafil
Output: CÓ

Ví dụ 2:
original_input: sex
identified_ingredient: Sildenafil
Output: KHÔNG

original_input: {original_input}
identified_ingredient: {identified_ingredient}
Output:
"""


PROMPT_GOC_RUT_GON = """
Bạn là một Dược sĩ lâm sàng AI chuyên nghiệp...
(Nội dung prompt này không đổi, giữ nguyên như cũ)
...
"""

# --- 2. CÁC HÀM XỬ LÝ (Cache) ---
@st.cache_resource
def get_model():
    print("--- Khởi tạo model AI ---")
    return genai.GenerativeModel('gemini-2.5-flash-lite')

@st.cache_data
def get_drug_info(drug_name):
    """
    Thực hiện quy trình tra cứu 3 BƯỚC (Nhận diện -> Xác thực -> Phân tích)
    và cache kết quả.
    """
    print(f"--- Bắt đầu quy trình tra cứu cho: {drug_name} ---")
    model = get_model()

    # --- BƯỚC 1A: NHẬN DIỆN SƠ BỘ ---
    prompt_nhan_dien_final = PROMPT_NHAN_DIEN.format(drug_name=drug_name)
    response_nhan_dien = model.generate_content(prompt_nhan_dien_final)
    hoat_chat_goc = response_nhan_dien.text.strip()

    # Nếu bước 1 trả về lỗi ngay, dừng lại
    if not hoat_chat_goc or "LỖI" in hoat_chat_goc:
        return f"❌ Lỗi: '{drug_name}' không được nhận dạng là một tên thuốc hợp lệ."

    # --- BƯỚC 1B: XÁC THỰC CHÉO ---
    print(f"--- Đã nhận diện sơ bộ: {hoat_chat_goc}. Bắt đầu xác thực... ---")
    prompt_xac_thuc_final = PROMPT_XAC_THUC.format(original_input=drug_name, identified_ingredient=hoat_chat_goc)
    response_xac_thuc = model.generate_content(prompt_xac_thuc_final)
    xac_thuc_text = response_xac_thuc.text.strip().upper()

    # --- BƯỚC 1C: QUYẾT ĐỊNH ---
    # Chỉ tiếp tục nếu câu trả lời xác thực là "CÓ"
    if "CÓ" not in xac_thuc_text:
        print(f"--- Xác thực thất bại cho '{drug_name}' và '{hoat_chat_goc}'. ---")
        return f"❌ Lỗi: '{drug_name}' không được nhận dạng là một tên thuốc hợp lệ."
    
    print(f"--- Xác thực thành công! Bắt đầu phân tích chi tiết. ---")
    # --- BƯỚC 2: PHÂN TÍCH CHI TIẾT ---
    full_prompt = f"{PROMPT_GOC_RUT_GON}\n\nHãy tra cứu và trình bày thông tin cho thuốc sau đây: **{hoat_chat_goc}**"
    response_phan_tich = model.generate_content(full_prompt)
    
    final_response = f"✅ Đã xác thực hoạt chất: **{hoat_chat_goc}**\n\n---\n\n{response_phan_tich.text}"
    return final_response


# --- 3. GIAO DIỆN VÀ LOGIC CHÍNH ---
st.title("Dược Điển AI")
# (Phần còn lại không thay đổi)
st.caption("Dự án được phát triển bởi group CÂCK và AI")

drug_name_input = st.text_input("Nhập tên thuốc (biệt dược hoặc hoạt chất):")
lookup_button = st.button("Tra cứu")

if lookup_button:
    if not drug_name_input:
        st.warning("Vui lòng nhập tên thuốc trước khi tra cứu.")
    else:
        try:
            with st.spinner("Dược sĩ AI đang làm việc, vui lòng chờ..."):
                final_result = get_drug_info(drug_name_input)
            st.markdown(final_result)
        
        except ga_ex.PermissionDenied as e:
            st.error("🚫 Lỗi Xác Thực: Google API Key của bạn không hợp lệ hoặc đã bị vô hiệu hóa...")
            st.exception(e)
        except ga_ex.ResourceExhausted as e:
            st.error("🚦 Đã đạt giới hạn: Bạn đã gửi quá nhiều yêu cầu...")
            st.exception(e)
        except ValueError as e:
            if "safety setting" in str(e):
                st.error("🔒 Nội dung bị chặn: Yêu cầu của bạn có thể đã vi phạm chính sách an toàn...")
                st.exception(e)
            else:
                st.error(f"Lỗi Dữ Liệu: Có vấn đề với dữ liệu đầu vào hoặc đầu ra...")
                st.exception(e)
        except ga_ex.GoogleAPICallError as e:
            st.error("🌐 Lỗi Kết Nối: Không thể kết nối đến máy chủ của Google AI...")
            st.exception(e)
        except Exception as e:
            st.error("💥 Lỗi không xác định: Một sự cố không mong muốn đã xảy ra.")
            st.exception(e)

