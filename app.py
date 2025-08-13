# --- Thư viện cần thiết ---
import streamlit as st
import requests 
import json
import time
from datetime import datetime, timedelta

# --- Cấu hình AI (nếu có) ---
try:
    api_key = st.secrets["GOOGLE_API_KEY"]
    genai.configure(api_key=api_key)
    gemini_model = genai.GenerativeModel('gemini-pro')
    is_api_configured = True
except (KeyError, AttributeError):
    is_api_configured = False
    gemini_model = None


# --- Hàm Bước 1: Lấy mã RxCUI từ tên thuốc ---
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
            name_data = name_response.json()
            standard_name = name_data.get('propConceptGroup', {}).get('propConcept', [{}])[0].get('propValue', drug_name)
            return standard_name, rxcui
        else:
            return None, None
    except Exception as e:
        st.error(f"Lỗi khi gọi API RxNorm: {e}")
        return None, None

# --- Hàm Bước 2: Lấy thông tin chi tiết từ openFDA ---
def get_fda_data(drug_name):
    base_url = "https://api.fda.gov/drug/label.json"
    params = {'search': f'openfda.generic_name:"{drug_name}"', 'limit': 1}
    try:
        response = requests.get(base_url, params=params)
        response.raise_for_status()
        data = response.json()
        if 'results' in data and len(data['results']) > 0:
            return data['results'][0]
        else:
            return None
    except Exception as e:
        st.error(f"Lỗi khi gọi API openFDA: {e}")
        return None

# --- Hàm MỚI - Bước 3A: Tìm kiếm nghiên cứu trên PubMed ---
def get_recent_studies_from_pubmed(drug_name, num_studies=3):
    """Tìm các ID bài báo mới nhất từ PubMed."""
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    # Tìm kiếm các bài báo trong 2 năm gần đây
    two_years_ago = (datetime.now() - timedelta(days=730)).strftime('%Y/%m/%d')
    search_term = f'("{drug_name}"[Title/Abstract]) AND ("{two_years_ago}"[Date - Publication] : "3000"[Date - Publication])'
    
    search_params = {
        'db': 'pubmed',
        'term': search_term,
        'retmode': 'json',
        'retmax': num_studies,
        'sort': 'relevance'
    }
    try:
        response = requests.get(base_url + "esearch.fcgi", params=search_params)
        response.raise_for_status()
        data = response.json()
        return data.get('esearchresult', {}).get('idlist', [])
    except Exception as e:
        st.warning(f"Không thể tìm kiếm trên PubMed: {e}")
        return []

# --- Hàm MỚI - Bước 3B: Lấy chi tiết và tóm tắt bằng Gemini ---
def summarize_studies_with_gemini(pmids):
    """Lấy chi tiết bài báo và dùng Gemini để tóm tắt."""
    if not pmids or not gemini_model:
        return "Không tìm thấy nghiên cứu mới hoặc API của Gemini chưa được cấu hình."

    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    summary_params = {'db': 'pubmed', 'id': ",".join(pmids), 'retmode': 'xml'}
    
    try:
        response = requests.get(base_url + "efetch.fcgi", params=summary_params)
        response.raise_for_status()
        # Do API PubMed trả về XML hơi phức tạp, ta chỉ lấy văn bản thô để AI xử lý
        raw_text = response.text

        prompt = f"""
        Dưới đây là dữ liệu thô từ một vài bài báo trên PubMed.
        Hãy đọc và tóm tắt những phát hiện chính từ các nghiên cứu này thành một vài gạch đầu dòng bằng Tiếng Việt.
        Tập trung vào kết quả, không cần mô tả phương pháp.
        
        Dữ liệu:
        {raw_text[:4000]}
        """
        
        ai_response = gemini_model.generate_content(prompt)
        return ai_response.text
    except Exception as e:
        return f"Lỗi khi tóm tắt bằng AI: {e}"

# --- Giao diện chính ---
st.set_page_config(page_title="Dược Điển AI - Bước 3", page_icon="💊", layout="wide")
st.title("💊 Dược Điển AI - Bước 3: Tích hợp AI tóm tắt")
st.write("Phát triển bởi group CÂCK và cộng sự AI.")

ten_thuoc = st.text_input("Nhập tên thuốc (tên gốc hoặc biệt dược):", placeholder="Ví dụ: Lipitor, Paracetamol...")

if st.button("Tra cứu thuốc"):
    if not ten_thuoc:
        st.warning("Vui lòng nhập tên thuốc.")
    else:
        # Bước 1
        standard_name, rxcui = get_rxcui_from_name(ten_thuoc)
        if not rxcui:
            st.error(f"Không tìm thấy thuốc '{ten_thuoc}'.")
        else:
            st.success(f"Bước 1: Tìm thấy tên gốc '{standard_name}' (RxCUI: {rxcui})")
            time.sleep(1)
            
            # Bước 2
            fda_data = get_fda_data(standard_name)
            if fda_data:
                st.success(f"Bước 2: Đã tìm thấy dữ liệu từ FDA!")
                time.sleep(1)

                # Bước 3
                with st.spinner("Bước 3: Đang tìm và tóm tắt các nghiên cứu mới nhất..."):
                    pmids = get_recent_studies_from_pubmed(standard_name)
                    summary = summarize_studies_with_gemini(pmids)
                st.success("Bước 3 hoàn thành!")

                # Hiển thị kết quả
                st.subheader("Kết quả tra cứu sơ bộ")
                with st.expander("Chỉ định (từ FDA)"):
                    st.markdown(fda_data.get('indications_and_usage', ['Không có.'])[0])
                with st.expander("Liều dùng (từ FDA)"):
                    st.markdown(fda_data.get('dosage_and_administration', ['Không có.'])[0])
                with st.expander("Nghiên cứu mới nhất (từ PubMed + AI)"):
                    st.markdown(summary)
            else:
                st.error(f"Không tìm thấy dữ liệu chi tiết cho '{standard_name}' trên openFDA.")
