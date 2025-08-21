import streamlit as st
from core.init_app import init_app
from core.state import init_session_state
from ui.main import render_lookup_page, render_prescription_analysis_page
from ui.sidebar import render_sidebar

st.set_page_config(page_title="Dược Điển AI", page_icon="💊", layout="wide")

# Khởi tạo Firebase và session state
init_app()
init_session_state()

# Sidebar
render_sidebar()

# Nội dung chính
st.title("Dược Điển AI")
if st.session_state.current_page == "Tra cứu Dược điển":
    render_lookup_page()
elif st.session_state.current_page == "Phân tích Đơn thuốc":
    render_prescription_analysis_page()
