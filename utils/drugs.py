import streamlit as st
from utils.prompts import PROMPT_NHAN_DIEN, PROMPT_REGULAR, PROMPT_SUMMARY
from utils.models import get_lookup_model, get_pro_model
from utils.pubmed import search_pubmed

def get_drug_info_from_api(drug_name, is_pro_user=False):
    identifier_model = get_lookup_model()
    prompt_nhan_dien_final = PROMPT_NHAN_DIEN.format(drug_name=drug_name)
    response_nhan_dien = identifier_model.generate_content(prompt_nhan_dien_final)
    response_text = response_nhan_dien.text

    try:
        hoat_chat_goc = response_text.split("Output:")[1].strip()
    except IndexError:
        hoat_chat_goc = response_text.strip()

    if hoat_chat_goc == "INVALID" or not hoat_chat_goc:
        return f"❌ Lỗi: '{drug_name}' không được nhận dạng.", None

    analysis_model = get_lookup_model()
    generation_config = {"max_output_tokens": 8192, "temperature": 0.6}
    full_prompt = f"{PROMPT_REGULAR}\n\nHãy tra cứu và trình bày thông tin cho thuốc sau đây: **{hoat_chat_goc}**"
    response_phan_tich = analysis_model.generate_content(full_prompt, generation_config=generation_config)
    base_response_text = response_phan_tich.text
    final_response = f"✅ Hoạt chất đã nhận diện: **{hoat_chat_goc}**\n\n---\n\n{base_response_text}"

    if is_pro_user:
        section_11_content = "\n\n---\n\n**11. Phân tích các Nghiên cứu Lâm sàng nổi bật (trong 2 năm gần đây):**\n"
        try:
            with st.spinner("Người dùng Pro: Đang truy vấn API của PubMed..."):
                search_context = search_pubmed(hoat_chat_goc)
                summary_prompt_final = PROMPT_SUMMARY.format(drug_name=hoat_chat_goc, search_results=search_context)
                summary_model = get_pro_model()
                summary_response = summary_model.generate_content(summary_prompt_final, generation_config=generation_config)
                section_11_content += summary_response.text
        except Exception as e:
            st.warning(f"Lỗi khi xử lý thông tin từ PubMed: {e}")
            section_11_content += "Đã xảy ra lỗi khi cố gắng tóm tắt dữ liệu từ PubMed."
        final_response += section_11_content

    return final_response, hoat_chat_goc


# ----------------------------
# Wrapper để page.py sử dụng
# ----------------------------
def lookup_drug(drug_name: str, is_pro_user: bool = False):
    """Wrapper để dùng trong UI"""
    result, _ = get_drug_info_from_api(drug_name, is_pro_user)
    return result