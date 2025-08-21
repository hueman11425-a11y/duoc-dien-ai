import streamlit as st
from utils import firebase

def handle_sidebar_click(drug_item):
    """
    Xử lý click thuốc từ lịch sử hoặc bộ sưu tập.
    drug_item có thể là string hoặc dict với key 'title'.
    """
    st.session_state.current_page = "Tra cứu Dược điển"
    
    # Lấy tên thuốc
    drug_name = drug_item.get("title") if isinstance(drug_item, dict) else str(drug_item)
    st.session_state.lookup_input_field = drug_name

    # Lấy kết quả đã lưu
    if st.session_state.user_info:
        result = firebase.load_user_result(st.session_state.firebase_db, st.session_state.user_info, drug_name)
    else:
        result = st.session_state.guest_cache.get(drug_name)
    
    st.session_state.query_result = result if result else "Không tìm thấy kết quả đã lưu."
