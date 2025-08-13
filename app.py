# --- Thư viện cần thiết ---
import streamlit as st
import google.generativeai as genai
import time
import json # <<< THÊM THƯ VIỆN JSON

# --- Cấu hình và khởi tạo mô hình AI ---
try:
    api_key = st.secrets["GOOGLE_API_KEY"]
    genai.configure(api_key=api_key)
    model = genai.GenerativeModel('gemini-2.5-pro')
    is_api_configured = True
except (KeyError, AttributeError):
    is_api_configured = False
except Exception as e:
    is_api_configured = False

# --- Cấu hình chung cho việc gọi AI ---
safety_settings = {
    "HARM_CATEGORY_HARASSMENT": "BLOCK_NONE",
    "HARM_CATEGORY_HATE_SPEECH": "BLOCK_NONE",
    "HARM_CATEGORY_SEXUALLY_EXPLICIT": "BLOCK_NONE",
    "HARM_CATEGORY_DANGEROUS_CONTENT": "BLOCK_NONE",
}
generation_config = {
    "max_output_tokens": 8192,
}

# --- Hàm gọi AI cho từng mục riêng lẻ ---
def get_drug_info_section(drug_name, section_name, section_prompt):
    """Hàm này yêu cầu AI trả về dữ liệu dưới dạng JSON."""
    
    # ===== KỸ THUẬT 3: YÊU CẦU ĐỊNH DẠNG JSON =====
    full_prompt = f"""
Nhiệm vụ: Trích xuất thông tin dược lý cho thuốc '{drug_name}' và trả về dưới dạng một đối tượng JSON.
Mục thông tin cần trích xuất: '{section_name}'.
Yêu cầu cho mục này: '{section_prompt}'.

QUAN TRỌNG: Chỉ trả về một đối tượng JSON hợp lệ duy nhất có cấu trúc sau:
{{
  "data": "Nội dung thông tin bằng Tiếng Việt ở đây"
}}
Không thêm bất kỳ văn bản, giải thích, hay định dạng markdown nào khác ngoài đối tượng JSON.
"""
    try:
        response = model.generate_content(
            full_prompt,
            generation_config=generation_config,
            safety_settings=safety_settings
        )
        
        # Xử lý kết quả trả về dạng JSON
        response_text = response.text.strip()
        
        # Đảm bảo văn bản trả về là một JSON hợp lệ
        # Đôi khi AI có thể trả về ```json\n{...}\n```, cần loại bỏ nó
        if response_text.startswith("```json"):
            response_text = response_text[7:-3].strip()

        data_obj = json.loads(response_text)
        return data_obj.get("data", "*Lỗi: AI trả về JSON nhưng thiếu key 'data'.*")

    except json.JSONDecodeError:
        return "*Lỗi: AI đã không trả về một định dạng JSON hợp lệ.*"
    except ValueError:
        return "*Lỗi: Phản hồi cho mục này đã bị chặn bởi bộ lọc an toàn.*"
    except Exception as e:
        return f"*Lỗi khi gọi AI: {e}*"

# --- Xây dựng giao diện ứng dụng với Streamlit (Không thay đổi) ---
st.set_page_config(page_title="Dược Điển AI", page_icon="💊", layout="wide")
st.title("💊 Dược Điển AI - Tra Cứu Dược Lý Thông Minh")
st.write("Cung cấp thông tin thuốc nhanh chóng, đáng tin cậy cho chuyên gia y tế. Phát triển bởi group CÂCK và cộng sự AI.")

sections = {
    "1. Tên thuốc": "Liệt kê tên gốc (in đậm) và các tên biệt dược phổ biến.",
    "2. Nhóm thuốc": "Nêu rõ phân loại dược lý.",
    "3. Cơ chế": "Giải thích rõ ràng, súc tích cách thuốc hoạt động.",
    "4. Dược động học (ADME)": "Trình bày đủ 4 mục con: Hấp thu (Absorption), Phân bố (Distribution), Chuyển hóa (Metabolism), Thải trừ (Excretion).",
    "5. Chỉ định": "Sử dụng danh sách gạch đầu dòng cho các chỉ định đã được cấp phép.",
    "6. Chống chỉ định": "Sử dụng danh sách gạch đầu dòng cho các trường hợp tuyệt đối không được dùng thuốc.",
    "7. Tương tác thuốc": "Liệt kê các tương tác quan trọng, giải thích ngắn gọn hậu quả.",
    "8. Tác dụng phụ": "Phân loại rõ ràng theo tần suất nếu có thể: Thường gặp, Ít gặp, Hiếm gặp.",
    "9. Lưu ý lâm sàng & Theo dõi": "Những cảnh báo quan trọng cho bác sĩ/dược sĩ và các xét nghiệm cần theo dõi khi dùng thuốc.",
    "10. Liều dùng": "Ghi rõ liều cho các chỉ định và đối tượng khác nhau nếu có thông tin.",
    "11. Nghiên cứu mới nhất": "Tóm tắt ngắn gọn 3-5 nghiên cứu nổi bật trong 1-2 năm gần đây, nêu rõ kết quả chính và tên tạp chí công bố."
}

if not is_api_configured:
    st.error("LỖI CẤU HÌNH: Google API Key chưa được thiết lập trong Streamlit Secrets. Vui lòng liên hệ quản trị viên.")
else:
    ten_thuoc = st.text_input("Nhập tên thuốc (tên gốc hoặc biệt dược):", placeholder="Ví dụ: Atorvastatin hoặc Lipitor")

    if st.button("Tra cứu Thông tin Thuốc"):
        if not ten_thuoc:
            st.warning("Vui lòng nhập tên thuốc cần tra cứu.")
        else:
            st.divider()
            st.subheader(f"Báo cáo chi tiết về {ten_thuoc}")
            
            for section_name, section_prompt in sections.items():
                with st.spinner(f"Đang lấy thông tin mục: {section_name}..."):
                    with st.expander(f"**{section_name}**", expanded=True):
                        result = get_drug_info_section(ten_thuoc, section_name, section_prompt)
                        st.markdown(result)
                
                time.sleep(3)
            
            st.divider()
            st.success("Hoàn tất tra cứu!")
            st.markdown("*Lưu ý: Thông tin trên chỉ mang tính chất tham khảo và không thể thay thế cho chẩn đoán, tư vấn và chỉ định của chuyên gia y tế. Luôn tham khảo ý kiến bác sĩ hoặc dược sĩ trước khi sử dụng bất kỳ loại thuốc nào.*")
