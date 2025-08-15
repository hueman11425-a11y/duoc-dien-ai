# Hàm tra cứu thuốc + PubMed
import streamlit as st
from datetime import date, timedelta
from Bio import Entrez
import time
from core.config import AppConfig

config = AppConfig()

@st.cache_data(ttl=3600)
def search_pubmed(drug_name: str) -> str:
    Entrez.email = "duocdien.ai.project@example.com"
    api_key = st.secrets.get("api_keys", {}).get("pubmed")
    if api_key:
        Entrez.api_key = api_key

    today = date.today()
    two_years_ago = today - timedelta(days=730)
    date_filter = f'AND ("{two_years_ago:%Y/%m/%d}"[Date - Publication] : "{today:%Y/%m/%d}"[Date - Publication])'

    search_term = f'"{drug_name}"[Title/Abstract] AND ("clinical trial"[Publication Type] OR "systematic review"[Publication Type]) {date_filter}'
    try:
        ids = Entrez.read(Entrez.esearch(db="pubmed", term=search_term, retmax="5", sort="relevance"))["IdList"]
        if not ids:
            return "Không tìm thấy bài báo phù hợp."

        records = Entrez.efetch(db="pubmed", id=ids, rettype="medline", retmode="text").read()
        context = ""
        for article in records.strip().split("\n\n"):
            title = _extract_field(article, "TI  - ")
            abstract = _extract_field(article, "AB  - ")
            journal = _extract_field(article, "JT  - ")
            pub_date = _extract_field(article, "DP  - ")[:4]
            pmid = _extract_field(article, "PMID- ")
            context += f"- Tiêu đề: {title}\n- Tạp chí: {journal}\n- Năm: {pub_date}\n- Tóm tắt: {abstract}\n- PMID: {pmid}\n\n"
            time.sleep(0.1)
        return context
    except Exception as e:
        return f"Lỗi PubMed: {e}"

def _extract_field(text: str, prefix: str) -> str:
    return next((line[len(prefix):] for line in text.split('\n') if line.startswith(prefix)), "N/A")

@st.cache_data(ttl="6h")
def get_drug_info(drug_name: str, is_pro_user: bool = False) -> str:
    model_id = config.models["regular"]
    hoat_chat_goc = _identify_active_ingredient(model_id, config.prompts["nhan_dien"], drug_name)
    if hoat_chat_goc == "INVALID":
        return f"❌ Lỗi: '{drug_name}' không được nhận dạng."

    analysis_model = config.models["pro"] if is_pro_user else model_id
    prompt = config.prompts["pro"] if is_pro_user else config.prompts["regular"]
    base_response = analysis_model.generate_content(
        f"{prompt}\n\nThông tin thuốc: **{hoat_chat_goc}**",
        generation_config={"max_output_tokens": 8192, "temperature": 0.6}
    ).text

    result = f"Hoạt chất: {hoat_chat_goc}\n{base_response}"
    if is_pro_user:
        pubmed_summary = search_pubmed(hoat_chat_goc)
        summary_response = config.models["pro"].generate_content(
            config.prompts["summary"].format(drug_name=hoat_chat_goc, search_results=pubmed_summary),
            generation_config={"max_output_tokens": 8192, "temperature": 0.6}
        ).text
        result += f"\n\n---\n\n**11. Nghiên cứu lâm sàng gần đây:**\n{summary_response}"
    return result

def _identify_active_ingredient(model, prompt_template, drug_name: str) -> str:
    response = model.generate_content(prompt_template.format(drug_name=drug_name))
    try:
        return response.text.split("Output:")[1].strip()
    except IndexError:
        return response.text.strip()
