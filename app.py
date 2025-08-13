import streamlit as st
import google.generativeai as genai
import google.generativeai.types as genai_types
import google.api_core.exceptions as ga_ex

# --- 1. KHỞI TẠO TRẠNG THÁI PHIÊN (SESSION STATE) ---
if 'history' not in st.session_state:
    st.session_state.history = []

# --- 2. CẤU HÌNH VÀ PROMPTS ---
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

# --- 3. CÁC HÀM XỬ LÝ (Cache) ---
@st.cache_resource
def get_model():
    return genai.GenerativeModel('gemini-2.5-flash-lite')

@st.cache_data
def get_drug_info(drug_name):
    model = get_model()
    prompt_nhan_dien_final = PROMPT_NHAN_DIEN.format(drug_name=drug_name)
    response_nhan_dien = model.generate_content(prompt_nhan_dien_final)
    hoat_chat_goc = response_nhan_dien.text.strip()
    if not hoat_chat_goc or "LỖI" in hoat_chat_goc:
        return f"❌ Lỗi: '{drug_name}' không được nhận dạng là một tên thuốc hợp lệ."
    prompt_xac_thuc_final = PROMPT_XAC_THUC.format(original_input=drug_name, identified_ingredient=hoat_chat_goc)
    response_xac_thuc = model.generate_content(prompt_xac_thuc_final)
    xac_thuc_text = response_xac_thuc.text.strip().upper()
    if "CÓ" not in xac_thuc_text:
        return f"❌ Lỗi: '{drug_name}' không được nhận dạng là một tên thuốc hợp lệ."
    full_prompt = f"{PROMPT_GOC_RUT_GON}\n\nHãy tra cứu và trình bày thông tin cho thuốc sau đây: **{hoat_chat_goc}**"
    response_phan_tich = model.generate_content(full_prompt)
    final_response = f"✅ Đã xác thực hoạt chất: **{hoat_chat_goc}**\n\n---\n\n{response_phan_tich.text}"
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
        st.error("🚫 Lỗi Xác Thực: Google API Key của bạn không hợp lệ hoặc đã bị vô hiệu hóa. Vui lòng kiểm tra lại trong file `.streamlit/secrets.toml`.")
        st.exception(e)
    except ga_ex.ResourceExhausted as e:
        st.error("🚦 Đã đạt giới hạn: Bạn đã gửi quá nhiều yêu cầu trong một thời gian ngắn. Vui lòng chờ vài phút rồi thử lại.")
        st.exception(e)
    except ValueError as e:
        if "safety setting" in str(e):
            st.error("🔒 Nội dung bị chặn: Yêu cầu của bạn có thể đã vi phạm chính sách an toàn của Google. Vui lòng thử lại với một tên thuốc khác.")
            st.exception(e)
        else:
            st.error(f"Lỗi Dữ Liệu: Có vấn đề với dữ liệu đầu vào hoặc đầu ra. Chi tiết: {e}")
            st.exception(e)
    except ga_ex.GoogleAPICallError as e:
        st.error("🌐 Lỗi Kết Nối: Không thể kết nối đến máy chủ của Google AI. Vui lòng kiểm tra lại kết nối mạng của bạn.")
        st.exception(e)
    except Exception as e:
        st.error("💥 Lỗi không xác định: Một sự cố không mong muốn đã xảy ra.")
        st.exception(e)

# --- 5. GIAO DIỆN VÀ LOGIC CHÍNH ---
st.title("Dược Điển AI Closed Beta 💊")
st.caption("Dự án được phát triển bởi group CÂCK và AI")

# --- HIỂN THỊ LỊCH SỬ TRÊN SIDEBAR ---
st.sidebar.header("Lịch sử tra cứu")
if not st.session_state.history:
    st.sidebar.info("Chưa có thuốc nào được tra cứu.")
else:
    for drug in st.session_state.history:
        if st.sidebar.button(drug, key=f"history_{drug}", use_container_width=True):
            run_lookup(drug)

# --- KHUNG GÓP Ý MỚI TRÊN SIDEBAR ---
with st.sidebar.container(border=True):
    st.write("**Bạn có ý tưởng để cải thiện ứng dụng?**")
    st.link_button(
        "Gửi phản hồi ngay!",
        url="https://forms.gle/M44GDS4hJ7LpY7b98", # <-- LINK CỦA BẠN ĐÂY
        help="Mở form góp ý trong một tab mới"
    )

# --- KHU VỰC NHẬP LIỆU CHÍNH ---
drug_name_input = st.text_input("Nhập tên thuốc (biệt dược hoặc hoạt chất):", key="main_input")
lookup_button = st.button("Tra cứu")

if lookup_button:
    if not drug_name_input:
        st.warning("Vui lòng nhập tên thuốc trước khi tra cứu.")
    else:
        run_lookup(drug_name_input)
