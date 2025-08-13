# --- Thư viện cần thiết ---
import streamlit as st
import requests # Thư viện để gọi API bên ngoài
import json

# --- Hàm mới: Lấy mã RxCUI từ tên thuốc ---
def get_rxcui_from_name(drug_name):
    """
    Hàm này kết nối đến API RxNorm của NIH để tìm mã định danh (RxCUI)
    và tên gốc chuẩn hóa từ một tên thuốc bất kỳ.
    """
    # URL của RxNorm API
    base_url = "https://rxnav.nlm.nih.gov/REST/rxcui.json"
    # Tạo tham số cho yêu cầu API
    params = {'name': drug_name, 'search': 1}
    
    try:
        # Gửi yêu cầu đến API
        response = requests.get(base_url, params=params)
        response.raise_for_status()  # Báo lỗi nếu có vấn đề về kết nối (lỗi 4xx hoặc 5xx)
        
        data = response.json()
        
        # Lấy mã RxCUI từ kết quả. Thường kết quả đầu tiên là chính xác nhất.
        id_group = data.get('idGroup', {})
        if id_group and 'rxnormId' in id_group:
            rxcui = id_group['rxnormId'][0]
            # Sau khi có RxCUI, ta có thể lấy tên chuẩn hóa
            name_url = f"https://rxnav.nlm.nih.gov/REST/rxcui/{rxcui}/property.json?propName=RxNorm%20Name"
            name_response = requests.get(name_url)
            name_response.raise_for_status()
            name_data = name_response.json()
            standard_name = name_data.get('propConceptGroup', {}).get('propConcept', [{}])[0].get('propValue', drug_name)
            
            return standard_name, rxcui
        else:
            return None, None # Không tìm thấy
            
    except requests.exceptions.RequestException as e:
        st.error(f"Lỗi kết nối đến RxNorm API: {e}")
        return None, None
    except json.JSONDecodeError:
        st.error("Lỗi xử lý dữ liệu từ RxNorm API.")
        return None, None

# --- Giao diện ứng dụng Streamlit ---
st.set_page_config(page_title="Dược Điển AI - Bước 1", page_icon="💊", layout="wide")
st.title("💊 Dược Điển AI - Bước 1: Kiểm tra 'Trạm phiên dịch'")
st.write("Phát triển bởi group CÂCK và cộng sự AI.")

ten_thuoc = st.text_input("Nhập tên thuốc (tên gốc hoặc biệt dược):", placeholder="Ví dụ: Lipitor, Paracetamol, Augmentin...")

if st.button("Kiểm tra Tên thuốc"):
    if not ten_thuoc:
        st.warning("Vui lòng nhập tên thuốc.")
    else:
        with st.spinner(f"Đang 'phiên dịch' tên thuốc '{ten_thuoc}'..."):
            standard_name, rxcui = get_rxcui_from_name(ten_thuoc)
            
            if rxcui:
                st.success("Tìm thấy thông tin chuẩn hóa!")
                st.markdown(f"**Tên bạn nhập:** `{ten_thuoc}`")
                st.markdown(f"**Tên gốc chuẩn:** `{standard_name}`")
                st.markdown(f"**Mã định danh (RxCUI):** `{rxcui}`")
                st.info("Tuyệt vời! 'Trạm phiên dịch' hoạt động tốt. Chúng ta đã sẵn sàng cho bước tiếp theo.")
            else:
                st.error(f"Không tìm thấy mã định danh cho thuốc '{ten_thuoc}'. Vui lòng kiểm tra lại tên thuốc.")
