# Cấu hình ứng dụng & load prompt
import streamlit as st
import google.generativeai as genai
from pathlib import Path

class AppConfig:
    """Cấu hình ứng dụng, gồm model, API keys, và prompt."""
    def __init__(self):
        try:
            genai.configure(api_key=st.secrets["GOOGLE_API_KEY"])
        except KeyError:
            st.error("Vui lòng cấu hình GOOGLE_API_KEY.")
            st.stop()

        self.prompts = {
            "nhan_dien": load_prompt("prompts/prompt_nhandien.txt"),
            "regular": load_prompt("prompts/prompt_regular.txt"),
            "pro": load_prompt("prompts/prompt_pro.txt"),
            "summary": load_prompt("prompts/prompt_summary.txt"),
        }

        self.models = {
            "regular": genai.GenerativeModel(st.secrets.get("models", {}).get("regular", "gemini-2.5-flash-lite")),
            "pro": genai.GenerativeModel(st.secrets.get("models", {}).get("pro", "gemini-pro")),
        }

def load_prompt(file_path: str) -> str:
    path = Path(file_path)
    if not path.exists():
        st.error(f"Không tìm thấy file prompt: {file_path}")
        st.stop()
    return path.read_text(encoding="utf-8")
