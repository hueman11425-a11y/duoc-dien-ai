import streamlit as st
from datetime import date, timedelta
from Bio import Entrez
import time

@st.cache_data(ttl=3600)
def search_pubmed(drug_name):
    Entrez.email = "duocdien.ai.project@example.com"
    api_key = st.secrets.get("api_keys", {}).get("pubmed")
    if api_key:
        Entrez.api_key = api_key

    today = date.today()
    two_years_ago = today - timedelta(days=730)
    date_filter = f'AND ("{two_years_ago.strftime("%Y/%m/%d")}"[Date - Publication] : "{today.strftime("%Y/%m/%d")}"[Date - Publication])'
    search_term = f'"{drug_name}"[Title/Abstract] AND ("clinical trial"[Publication Type] OR "systematic review"[Publication Type]) {date_filter}'

    try:
        handle = Entrez.esearch(db="pubmed", term=search_term, retmax="5", sort="relevance")
        record = Entrez.read(handle)
        handle.close()
        id_list = record["IdList"]
        if not id_list:
            return "Không tìm thấy bài báo phù hợp nào trong 2 năm gần đây trên PubMed."
        
        handle = Entrez.efetch(db="pubmed", id=id_list, rettype="medline", retmode="text")
        records_text = handle.read()
        handle.close()

        context = ""
        articles = records_text.strip().split("\n\n")
        for article_text in articles:
            title = next((line[6:] for line in article_text.split('\n') if line.startswith("TI  - ")), "N/A")
            abstract = next((line[6:] for line in article_text.split('\n') if line.startswith("AB  - ")), "N/A")
            journal = next((line[6:] for line in article_text.split('\n') if line.startswith("JT  - ")), "N/A")
            pub_date = next((line[6:] for line in article_text.split('\n') if line.startswith("DP  - ")), "N/A")
            pmid = next((line[6:] for line in article_text.split('\n') if line.startswith("PMID- ")), "N/A")

            context += f"- Tiêu đề: {title}\n- Tạp chí: {journal}\n- Năm: {pub_date[:4]}\n- Tóm tắt: {abstract}\n- PMID: {pmid.strip()}\n\n"
            time.sleep(0.1)

        return context
    except Exception as e:
        return f"Đã xảy ra lỗi khi truy vấn API của PubMed: {e}"
