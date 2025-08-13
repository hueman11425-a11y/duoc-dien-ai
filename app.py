# --- Thư viện cần thiết ---
import streamlit as st
import google.generativeai as genai

# --- Cấu hình và khởi tạo mô hình AI ---
try:
    api_key = st.secrets["GOOGLE_API_KEY"]
    genai.configure(api_key=api_key)
    # ===== THAY ĐỔI 1: Nâng cấp model lên gemini-1.5-pro =====
    model = genai.GenerativeModel('gemini-2.5-pro')
    is_api_configured = True
except (KeyError, AttributeError):
    is_api_configured = False
except Exception as e:
    is_api_configured = False


# --- "Prompt Gốc" chính thức của dự án ---
PROMPT_GOC = """
Bạn là một Dược sĩ lâm sàng AI chuyên nghiệp và là chuyên gia trong việc tổng hợp thông tin y khoa. Nhiệm vụ của bạn là tra cứu và phân tích thông tin về một loại thuốc mà tôi cung cấp, sau đó trình bày kết quả **bằng Tiếng Việt**.

Hãy sử dụng toàn bộ kiến thức đã được huấn luyện của bạn từ các nguồn dữ liệu y khoa uy tín trên thế giới như sách giáo khoa (Goodman & Gilman's, Katzung's), các cơ sở dữ liệu mở (openFDA, WHO), và các tạp chí khoa học hàng đầu (PubMed, The Lancet, NEJM, Nature Reviews Drug Discovery, Pharmacological Reviews).

Khi tôi đưa tên một loại thuốc, bạn **PHẢI** trình bày kết quả theo đúng cấu trúc 11 mục sau đây, **sử dụng định dạng Markdown** để làm rõ các tiêu đề và danh sách:

**1. Tên thuốc:**
*(Liệt kê tên gốc (in đậm) và các tên biệt dược phổ biến)*

**2. Nhóm thuốc:**
*(Nêu rõ phân loại dược lý)*

**3. Cơ chế:**
*(Giải thích rõ ràng, súc tích cách thuốc hoạt động)*

**4. Dược động học (ADME):**
*(Trình bày dưới dạng danh sách 4 mục con)*
* **Hấp thu (Absorption):**
* **Phân bố (Distribution):**
* **Chuyển hóa (Metabolism):**
* **Thải trừ (Excretion):**

**5. Chỉ định:**
*(Sử dụng danh sách gạch đầu dòng cho các chỉ định đã được cấp phép)*

**6. Chống chỉ định:**
*(Sử dụng danh sách gạch đầu dòng cho các trường hợp tuyệt đối không được dùng thuốc)*

**7. Tương tác thuốc:**
*(Liệt kê các tương tác quan trọng, giải thích ngắn gọn hậu quả)*

**8. Tác dụng phụ:**
*(Phân loại rõ ràng theo tần suất nếu có thể: Thường gặp, Ít gặp, Hiếm gặp)*

**9. Lưu ý lâm sàng & Theo dõi:**
*(Những cảnh báo quan trọng cho bác sĩ/dược sĩ và các xét nghiệm cần theo dõi khi dùng thuốc)*

**10. Liều dùng:**
*(Ghi rõ liều cho các chỉ định và đối tượng khác nhau nếu có thông tin)*

**11. Nghiên cứu mới nhất:**
*(Tóm tắt ngắn gọn 3-5 nghiên cứu nổi bật trong 1-2 năm gần đây. Mỗi nghiên cứu trình bày trong 2-3 câu, nêu rõ logic khoa học, kết quả chính và tên tạp chí công bố)*

---
**QUY TẮC BẮT BUỘC:**
* Tuyệt đối KHÔNG được bịa đặt hay suy diễn thông tin.
* Nếu không tìm thấy dữ liệu cho mục nào, hãy ghi rõ: `Không có đủ dữ liệu đáng tin cậy.`
* Luôn ưu tiên thông tin được chấp thuận bởi FDA.
* **Quan trọng:** Ở cuối mỗi kết quả tra cứu, **luôn thêm câu sau**: *Lưu ý: Thông tin trên chỉ mang tính chất tham khảo và không thể thay thế cho chẩn đoán, tư vấn và chỉ định của chuyên gia y tế. Luôn tham khảo ý kiến bác sĩ hoặc dược sĩ trước khi sử dụng bất kỳ loại thuốc nào.*
"""

# --- Xây dựng giao diện ứng dụng với Streamlit ---
st.set_page_config(page_title="Dược Điển AI", page_icon="💊", layout="wide")

st.title("💊 Dược Điển AI - Tra Cứu Dược Lý Thông Minh")
# ===== THAY ĐỔI 2: Cập nhật tên nhóm phát triển =====
st.write("Cung cấp thông tin thuốc nhanh chóng, đáng tin cậy cho chuyên gia y tế. Phát triển bởi group CÂCK và cộng sự AI.")

if not is_api_configured:
    st.error("LỖI CẤU HÌNH: Google API Key chưa được thiết lập trong Streamlit Secrets. Vui lòng liên hệ quản trị viên.")
else:
    ten_thuoc = st.text_input("Nhập tên thuốc (tên gốc hoặc biệt dược):", placeholder="Ví dụ: Atorvastatin hoặc Lipitor")

    if st.button("Tra cứu Thông tin Thuốc"):
        if not ten_thuoc:
            st.warning("Vui lòng nhập tên thuốc cần tra cứu.")
        else:
            with st.spinner(f"Đang tổng hợp thông tin cho **{ten_thuoc}**... Quá trình này có thể mất vài chục giây."):
                try:
                    # Thiết lập cài đặt an toàn để cho phép các nội dung y khoa
                    safety_settings = {
                        "HARM_CATEGORY_HARASSMENT": "BLOCK_NONE",
                        "HARM_CATEGORY_HATE_SPEECH": "BLOCK_NONE",
                        "HARM_CATEGORY_SEXUALLY_EXPLICIT": "BLOCK_NONE",
                        "HARM_CATEGORY_DANGEROUS_CONTENT": "BLOCK_NONE",
                    }
                    
                    full_prompt = PROMPT_GOC + "\n\n" + f"Hãy tra cứu thông tin về thuốc sau: **{ten_thuoc}**"
                    
                    # Gọi API của Gemini với cài đặt an toàn
                    response = model.generate_content(
                        full_prompt,
                        safety_settings=safety_settings
                    )
                    
                    st.divider()
                    st.subheader(f"Báo cáo chi tiết về {ten_thuoc}")
                    st.markdown(response.text)
                    st.divider()

                except ValueError:
                     # Xử lý lỗi do bộ lọc an toàn khi không có response.text
                     st.error("Lỗi: Phản hồi từ AI đã bị chặn bởi bộ lọc an toàn. Điều này có thể xảy ra với các loại thuốc có thông tin nhạy cảm. Chúng tôi đang làm việc để cải thiện vấn đề này.")
                except Exception as e:
                    st.error(f"Đã có lỗi xảy ra trong quá trình gọi AI: {e}")

