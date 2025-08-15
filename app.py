import streamlit as st
from core.access_control import verify_code
from core.drug_lookup import get_drug_info
from ui.sidebar import render_sidebar
from ui.main import render_main_interface

# --- KIỂM TRA BẢO TRÌ ---
if st.secrets.get("maintenance_mode", False):
    st.set_page_config(page_title="Bảo trì", page_icon="🛠️")
    st.title("🛠️ Dược Điển AI đang được bảo trì")
    st.info(st.secrets.get("maintenance_message", "Ứng dụng đang được cập nhật."))
    st.stop()

# --- KHỞI TẠO PHIÊN ---
st.session_state.setdefault("history", [])
st.session_state.setdefault("pro_access", False)

# --- GIAO DIỆN ---
st.set_page_config(page_title="Dược Điển AI", page_icon="💊")
st.title("Dược Điển AI 💊")
st.caption("Dự án được phát triển bởi group Kiệt, Hào và AI")

render_sidebar(verify_code)
render_main_interface(get_drug_info)
