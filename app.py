# --- Thư viện cần thiết ---
import streamlit as st
import google.generativeai as genai
import requests 
import json
import time
from datetime import datetime, timedelta
import xml.etree.ElementTree as ET # <<< THÊM THƯ VIỆN XỬ LÝ XML

# --- Cấu hình AI (nếu có) ---
try:
    api_key = st.secrets["GOOGLE_API_KEY"]
    genai.configure(api_key=api_key)
    gemini_model = genai.GenerativeModel('gemini-1.5-pro')
    is_api_configured = True
except (KeyError, AttributeError):
    is_api_configured = False
    gemini_model = None
except Exception as e:
    is_api_configured = False
    gemini_model = None
    st.error(f"Lỗi khởi tạo AI model: {e}")

# --- Các hàm chức năng ---
# (Các hàm get_rxcui_from_name và get_fda_data không đổi)
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

def get_recent_studies_from_pubmed(drug_name, num_studies=3):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    two_years_ago = (datetime.now() - timedelta(days=730)).strftime('%Y/%m/%d')
    search_term = f'("{drug_name}"[Title/Abstract]) AND ("{two_years_ago}"[Date - Publication] : "3000"[Date - Publication])'
    search_params = {'db': 'pubmed', 'term': search_term, 'retmode': 'json', 'retmax': num_studies, 'sort': 'relevance'}
    try:
        response = requests.get(base_url + "esearch.fcgi", params=search_params)
        response.raise_for_status()
        data = response.json()
        return data.get('esearchresult', {}).get('idlist', [])
    except Exception as e:
        st.warning(f"Không thể tìm kiếm trên PubMed: {e}")
        return []

# ===== HÀM ĐƯỢC NÂNG CẤP: Xử lý XML trước khi tóm tắt =====
def summarize_studies_with_gemini(pmids):
    if not pmids:
        return "Không tìm thấy nghiên cứu mới trên PubMed."
    if not is_api_configured or not gemini_model:
        return "API của Gemini chưa được cấu hình, không thể tóm tắt."

    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    summary_params = {'db': 'pubmed', 'id': ",".join(pmids), 'retmode': 'xml'}
    
    try:
        response = requests.get(base_url + "efetch.fcgi", params=summary_params)
        response.raise_for_status()
        
        # --- BƯỚC DỌN DẸP DỮ LIỆU ---
        clean_text_parts = []
        root = ET.fromstring(response.content)
        for article in root.findall('.//PubmedArticle'):
            title_element = article.find('.//ArticleTitle')
            title = title_element.text if title_element is not None else "Không có tiêu đề"
            
            abstract_element = article.find('.//AbstractText')
            abstract = abstract_element.text if abstract_element is not None else "Không có tóm tắt."
            
            clean_text_parts.append(f"Tiêu đề: {title}\nTóm tắt: {abstract}\n---")
        
        clean_text = "\n".join(clean_text_parts)

        if not clean_text:
            return "Không thể trích xuất nội dung từ các bài báo PubMed."

        prompt = f"""
        Dưới đây là nội dung đã được làm sạch từ một vài bài báo trên PubMed.
        Hãy đọc và tóm tắt những phát hiện chính từ các nghiên cứu này thành một vài gạch đầu dòng bằng Tiếng Việt.
        Tập trung vào kết quả, không cần mô tả phương pháp.
        
        Nội dung:
        {clean_text}
        """
        
        ai_response = gemini_model.generate_content(prompt)
        return ai_response.text
    except Exception as e:
        return f"Lỗi khi tóm tắt bằng AI: {e}"

# --- Giao diện chính (không đổi) ---
st.set_page_config(page_title="Dược Điển AI - Bước 3", page_icon="💊", layout="wide")
st.title("💊 Dược Điển AI - Bước 3: Tích hợp AI tóm tắt")
st.write("Phát triển bởi group CÂCK và cộng sự AI.")

ten_thuoc = st.text_input("Nhập tên thuốc (tên gốc hoặc biệt dược):", placeholder="Ví dụ: Lipitor, Paracetamol...")

if st.button("Tra cứu thuốc"):
    if not ten_thuoc:
        st.warning("Vui lòng nhập tên thuốc.")
    else:
        standard_name, rxcui = get_rxcui_from_name(ten_thuoc)
        if not rxcui:
            st.error(f"Không tìm thấy thuốc '{ten_thuoc}'.")
        else:
            st.success(f"Bước 1: Tìm thấy tên gốc '{standard_name}' (RxCUI: {rxcui})")
            time.sleep(1)
            
            fda_data = get_fda_data(standard_name)
            if fda_data:
                st.success(f"Bước 2: Đã tìm thấy dữ liệu từ FDA!")
                time.sleep(1)

                with st.spinner("Bước 3: Đang tìm và tóm tắt các nghiên cứu mới nhất..."):
                    pmids = get_recent_studies_from_pubmed(standard_name)
                    summary = summarize_studies_with_gemini(pmids)
                st.success("Bước 3 hoàn thành!")

                st.subheader("Kết quả tra cứu sơ bộ")
                with st.expander("Chỉ định (từ FDA)"):
                    st.markdown(fda_data.get('indications_and_usage', ['Không có.'])[0])
                with st.expander("Liều dùng (từ FDA)"):
                    st.markdown(fda_data.get('dosage_and_administration', ['Không có.'])[0])
                with st.expander("Nghiên cứu mới nhất (từ PubMed + AI)"):
                    st.markdown(summary)
            else:
                st.error(f"Không tìm thấy dữ liệu chi tiết cho '{standard_name}' trên openFDA.")
