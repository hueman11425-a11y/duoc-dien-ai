# Khu vực nội dung chính
import streamlit as st

def render_main_interface(get_drug_info):
    drug_name = st.text_input("Nhập tên thuốc:")
    if st.button("Tra cứu"):
        if not drug_name:
            st.warning("Vui lòng nhập tên thuốc.")
        else:
            is_pro = st.session_state.get("pro_access", False)
            result = get_drug_info(drug_name, is_pro_user=is_pro)
            if not result.startswith("❌ Lỗi:"):
                st.markdown(result)
                _update_history(drug_name)
            else:
                st.error(result)

def _update_history(drug_name: str):
    if drug_name not in st.session_state.history:
        st.session_state.history.insert(0, drug_name)
        if len(st.session_state.history) > 10:
            st.session_state.history.pop()
