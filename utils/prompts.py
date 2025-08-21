import streamlit as st

def load_prompt(file_path):
    try:
        with open(file_path, "r", encoding="utf-8") as f:
            return f.read()
    except FileNotFoundError:
        st.error(f"LỖI: Không tìm thấy file prompt tại '{file_path}'.")
        st.stop()

PROMPT_NHAN_DIEN = load_prompt("prompts/nhandien.txt")
PROMPT_REGULAR = load_prompt("prompts/regular.txt")
PROMPT_PRO = load_prompt("prompts/pro.txt")
PROMPT_SUMMARY = load_prompt("prompts/summary.txt")
PROMPT_PRESCRIPTION = load_prompt("prompts/prescription.txt")
