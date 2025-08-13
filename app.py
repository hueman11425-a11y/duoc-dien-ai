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
    
PROMPT_NHAN_DIEN = """Nhiệm vụ của bạn là nhận diện (các) hoạt chất gốc từ tên thuốc được cung cấp.
- Nếu Input là một biệt dược đã biết, hãy trả về (các) hoạt chất gốc của nó.
- Nếu Input đã là một hoạt chất gốc, hãy trả về chính nó.
- Nếu Input KHÔNG PHẢI là tên thuốc (ví dụ: một từ thông thường, một câu vô nghĩa), bạn BẮT BUỘC phải trả về duy nhất từ "INVALID".
- Chỉ trả về tên hoạt chất hoặc từ "INVALID". Tuyệt đối không giải thích.

Ví dụ:
Input: Lipitor
Output: Atorvastatin
Input: Paracetamol
Output: Paracetamol
Input: Troysar AM
Output: Losartan, Amlodipine
Input: a cat
Output: INVALID
Input: sex
Output: INVALID

Input: {drug_name}
Output:
"""

PROMPT_GOC_RUT_GON = """
Bạn là một Dược sĩ lâm sàng AI chuyên nghiệp và là chuyên gia trong việc tổng hợp thông tin y khoa.
Nhiệm vụ của bạn là tra cứu và phân tích thông tin về một loại thuốc mà tôi cung cấp.
Hãy sử dụng toàn bộ kiến thức đã được huấn luyện của bạn từ các nguồn dữ liệu y khoa uy tín trên thế giới như sách giáo khoa (Goodman & Gilman's, Katzung's), các cơ sở dữ liệu mở (openFDA, WHO), và các tạp chí khoa học hàng đầu (PubMed, The Lancet, NEJM).

---
**LƯU Ý ĐẶC BIỆT KHI PHÂN TÍCH:**
- Nếu tên thuốc đầu vào chỉ có MỘT hoạt chất, hãy phân tích bình thường.
- Nếu tên thuốc đầu vào chứa NHIỀU hoạt chất (ví dụ: 'Losartan, Amlodipine'), hãy phân tích chúng như một **liệu pháp phối hợp**. Trong mỗi mục, hãy làm rõ vai trò của từng thành phần và cách chúng tác động qua lại nếu có. Đừng hỏi lại, hãy tiến hành phân tích ngay.
---

Khi tôi đưa tên một loại thuốc (luôn là tên gốc/hoạt chất), bạn PHẢI trình bày kết quả theo đúng cấu trúc 10 mục sau đây, sử dụng ngôn ngữ chuyên môn, chính xác và rõ ràng:

1.  **Tên thuốc:** (Tên gốc và các tên biệt dược phổ biến của sự phối hợp này)
2.  **Nhóm thuốc:** (Nêu nhóm của từng hoạt chất)
3.  **Cơ chế:** (Mô tả cơ chế của từng hoạt chất và lợi ích của việc phối hợp chúng)
4.  **Dược động học (ADME):** (Nêu các thông số chính cho cả hai hoạt chất nếu có khác biệt đáng kể)
5.  **Chỉ định:** (Các chỉ định được phê duyệt cho liệu pháp phối hợp này)
6.  **Chống chỉ định:**
7.  **Tương tác thuốc:**
8.  **Tác dụng phụ:** (Nêu các tác dụng phụ chung và riêng của từng thành phần)
9.  **Lưu ý lâm sàng & Theo dõi:**
10. **Liều dùng:** (Trình bày liều dùng cụ thể cho các chỉ định chính của thuốc phối hợp)

**QUY TẮC BẮT BUỘC:**
- Tuyệt đối KHÔNG được bịa đặt hay suy diễn thông tin.
- Nếu không tìm thấy dữ liệu cho mục nào, hãy ghi rõ: "Không có đủ dữ liệu đáng tin cậy."
- Luôn ưu tiên thông tin được chấp thuận bởi FDA.
"""

# --- 3. CÁC HÀM XỬ LÝ (Cache) ---
@st.cache_resource
def get_model():
    return genai.GenerativeModel('gemini-2.5-flash-lite')

@st.cache_data(ttl="6h") # Thêm ttl để dữ liệu tự làm mới sau 6 tiếng
def get_drug_info(drug_name):
    model = get_model()
    prompt_nhan_dien_final = PROMPT_NHAN_DIEN.format(drug_name=drug_name)
    response_nhan_dien = model.generate_content(prompt_nhan_dien_final)
    hoat_chat_goc = response_nhan_dien.text.strip()
    if hoat_chat_goc == "INVALID":
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
