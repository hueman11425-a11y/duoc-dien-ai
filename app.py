# --- Thư viện cần thiết ---
import streamlit as st
import requests 
import json

# --- Hàm Bước 1: Lấy mã RxCUI từ tên thuốc (đã có từ trước) ---
def get_rxcui_from_name(drug_name):
    base_url = "https://rxnav.nlm.nih.gov/REST/rxcui.json"
    params = {'name': drug_name, 'search': 1}
    try:
        response = requests.get(base_url, params=params)
        response.raise_for_status()
        data = response.json()
        id_group = data.get('idGroup', {})
        if id_group and 'rxnormId' in id_group:
            rxcui = id_group['rxnormId'][0]
            name_url = f"https://rxnav.nlm.nih.gov/REST/rxcui/{rxcui}/property.json?propName=RxNorm%20Name"
            name_response = requests.get(name_url)
            name_response.raise_for_status()
            name_data = name_response.json()
            standard_name = name_data.get('propConceptGroup', {}).get('propConcept', [{}])[0].get('propValue', drug_name)
            return standard_name, rxcui
        else:
            return None, None
    except requests.exceptions.RequestException as e:
        st.error(f"Lỗi kết nối đến RxNorm API: {e}")
        return None, None
    except json.JSONDecodeError:
        st.error("Lỗi xử lý dữ liệu từ RxNorm API.")
        return None, None

# --- Hàm MỚI - Bước 2: Lấy thông tin chi tiết từ openFDA ---
def get_fda_data(drug_name):
    """
    Hàm này kết nối đến API openFDA để tìm thông tin chi tiết về thuốc
    dựa trên tên gốc (generic name).
    """
    # URL của openFDA API cho nhãn thuốc
    base_url = "https://api.fda.gov/drug/label.json"
    # Tìm kiếm dựa trên tên gốc trong trường openfda của nhãn
    params = {'search': f'openfda.generic_name:"{drug_name}"', 'limit': 1}
    
    try:
        response = requests.get(base_url, params=params)
        response.raise_for_status()
        data = response.json()
        
        # Lấy kết quả đầu tiên
        if 'results' in data and len(data['results']) > 0:
            drug_info = data['results'][0]
            # Trả về toàn bộ thông tin của nhãn thuốc đầu tiên tìm thấy
            return drug_info
        else:
            return None # Không tìm thấy
            
    except requests.exceptions.RequestException as e:
        st.error(f"Lỗi kết nối đến openFDA API: {e}")
        return None
    except json.JSONDecodeError:
        st.error("Lỗi xử lý dữ liệu từ openFDA API.")
        return None

# --- Giao diện ứng dụng Streamlit ---
st.set_page_config(page_title="Dược Điển AI - Bước 2", page_icon="💊", layout="wide")
st.title("💊 Dược Điển AI - Bước 2: Kết nối 'Kho dữ liệu'")
st.write("Phát triển bởi group CÂCK và cộng sự AI.")

ten_thuoc = st.text_input("Nhập tên thuốc (tên gốc hoặc biệt dược):", placeholder="Ví dụ: Lipitor, Paracetamol, Augmentin...")

if st.button("Tra cứu thuốc"):
    if not ten_thuoc:
        st.warning("Vui lòng nhập tên thuốc.")
    else:
        # --- Thực hiện Bước 1 ---
        with st.spinner(f"Bước 1: Đang 'phiên dịch' tên thuốc '{ten_thuoc}'..."):
            standard_name, rxcui = get_rxcui_from_name(ten_thuoc)
        
        if not rxcui:
            st.error(f"Không tìm thấy mã định danh cho thuốc '{ten_thuoc}'. Vui lòng kiểm tra lại tên thuốc.")
        else:
            st.success(f"Bước 1 thành công: Tìm thấy tên gốc là '{standard_name}' (RxCUI: {rxcui})")
            
            # --- Thực hiện Bước 2 ---
            with st.spinner(f"Bước 2: Đang tìm kiếm '{standard_name}' trong kho dữ liệu của FDA..."):
                fda_data = get_fda_data(standard_name)

            if fda_data:
                st.success("Bước 2 thành công: Đã tìm thấy dữ liệu từ FDA!")
                
                # Hiển thị thử một vài mục để kiểm tra
                st.subheader("Dữ liệu thô từ FDA (để kiểm tra)")

                # Dùng st.expander để không làm rối giao diện
                with st.expander("Chỉ định (Indications and Usage)"):
                    # Dữ liệu trả về là một danh sách, ta lấy phần tử đầu tiên
                    indications = fda_data.get('indications_and_usage', ['Không có thông tin.'])[0]
                    st.markdown(indications)

                with st.expander("Chống chỉ định (Contraindications)"):
                    contraindications = fda_data.get('contraindications', ['Không có thông tin.'])[0]
                    st.markdown(contraindications)
                
                with st.expander("Liều dùng và Cách dùng (Dosage and Administration)"):
                    dosage = fda_data.get('dosage_and_administration', ['Không có thông tin.'])[0]
                    st.markdown(dosage)
                
            else:
                st.error(f"Không tìm thấy thông tin chi tiết cho '{standard_name}' trên openFDA.")
