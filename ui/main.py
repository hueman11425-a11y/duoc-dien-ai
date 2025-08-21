import streamlit as st
from utils import drugs, firebase, constants, prescription

def render_lookup_page():
    st.header("Tra cứu Dược điển 💊")
    is_logged_in = st.session_state.user_info is not None
    firebase_db = st.session_state.firebase_db
    user_info = st.session_state.user_info
    
    def run_lookup(drug_name):
        st.session_state.query_result = None
        is_pro = st.session_state.get("pro_access", False)
        with st.spinner(f"Đang tra cứu '{drug_name}'..."):
            api_result, identified_name = drugs.get_drug_info_from_api(drug_name, is_pro)
            st.session_state.query_result = api_result
            if identified_name and "không được nhận dạng" not in api_result:
                if is_logged_in:
                    updated_history = firebase.save_new_result(firebase_db, user_info, identified_name, api_result)
                    if updated_history is not None:
                        st.session_state.history = updated_history
                else:
                    st.session_state.guest_cache[identified_name] = api_result
                    if identified_name not in st.session_state.history:
                        st.session_state.history.insert(0, identified_name)
                        if len(st.session_state.history) > 10:
                            drug_to_remove = st.session_state.history.pop()
                            st.session_state.guest_cache.pop(drug_to_remove, None)

    search_query = st.text_input("Nhập tên thuốc:", key="lookup_input_field")
    
    if st.button("Tra cứu"):
        if search_query:
            run_lookup(search_query)
        else:
            st.warning("Vui lòng nhập tên thuốc trước khi tra cứu.")

    if st.session_state.query_result:
        result_to_display = st.session_state.query_result
        if result_to_display.startswith("❌ Lỗi:"):
            st.error(result_to_display)
        else:
            st.markdown(result_to_display)
            if is_logged_in:
                identified_name = result_to_display.split("**")[1] if "**" in result_to_display else None
                if identified_name:
                    st.markdown("---")
                    st.subheader(f"Lưu '{identified_name}' vào bộ sưu tập")
                    collections = st.session_state.get("collections", {})
                    if not collections:
                        st.info("Bạn chưa có bộ sưu tập nào.")
                    else:
                        col1, col2 = st.columns([2,1])
                        with col1:
                            sel_coll = st.selectbox("Chọn bộ sưu tập:", options=list(collections.keys()))
                        with col2:
                            st.write(""); st.write("")
                            if st.button("Thêm thuốc", use_container_width=True):
                                msg = firebase.add_drug_to_collection(firebase_db, user_info, sel_coll, identified_name)
                                st.toast(msg)
                                _, collections_new = firebase.load_user_data(firebase_db, user_info)
                                st.session_state.collections = collections_new
                                st.rerun()


def render_prescription_analysis_page():
    st.header("Phân tích Đơn thuốc 🩺")
    is_logged_in = st.session_state.user_info is not None
    
    if not is_logged_in:
        st.warning("Vui lòng đăng nhập để sử dụng tính năng này.")
        return
    
    st.subheader("Phần 1: Bối cảnh Bệnh nhân (Không bắt buộc)")
    col1, col2 = st.columns(2)
    with col1:
        conditions = st.text_area("Tình trạng bệnh lý nền:", height=150)
    with col2:
        allergies = st.text_area("Dị ứng thuốc đã biết:", height=150)
        
    st.subheader("Phần 2: Thông tin Đơn thuốc (Bắt buộc)")
    prescription_text = st.text_area("Dán nội dung đơn thuốc:", height=250)
    
    if st.button("Phân tích Đơn thuốc", type="primary"):
        if not prescription_text.strip():
            st.error("Vui lòng nhập thông tin đơn thuốc.")
        else:
            with st.spinner("AI đang phân tích..."):
                patient_context = f"- Bệnh lý nền: {conditions or 'Không có'}\n- Dị ứng: {allergies or 'Không có'}"
                result = prescription.get_prescription_analysis(
                    st.session_state.firebase_db, st.session_state.user_info,
                    patient_context, prescription_text
                )
                st.session_state.analysis_result = result
    
    if st.session_state.analysis_result:
        st.markdown("---")
        st.subheader("Kết quả Phân tích")
        if st.session_state.analysis_result.startswith("❌"):
             st.error(st.session_state.analysis_result)
        else:
            st.markdown(st.session_state.analysis_result)
