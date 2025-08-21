import streamlit as st
import auth
from utils import firebase, constants
from core.callbacks import handle_sidebar_click
from uuid import uuid4

def render_sidebar():
    with st.sidebar:
        # --- Form đăng nhập ---
        is_logged_in = auth.display_auth_forms(st.session_state.firebase_auth, st.session_state.firebase_db)

        # --- Load dữ liệu user ---
        if is_logged_in and not st.session_state.user_data_loaded:
            with st.spinner("Đang tải dữ liệu..."):
                history, collections = firebase.load_user_data(st.session_state.firebase_db, st.session_state.user_info)
                st.session_state.history = history
                st.session_state.collections = collections
            st.session_state.user_data_loaded = True

        st.markdown("---")
        st.header("Tính năng")
        page_options = ["Tra cứu Dược điển", "Phân tích Đơn thuốc"]

        def on_page_change():
            st.session_state.current_page = st.session_state.feature_radio
            st.session_state.query_result = None
            st.session_state.analysis_result = None

        current_page_index = page_options.index(st.session_state.current_page) if st.session_state.current_page in page_options else 0
        st.radio(
            "Chọn tính năng:", page_options,
            index=current_page_index,
            on_change=on_page_change,
            key='feature_radio',
            label_visibility="collapsed"
        )

        # --- Lịch sử ---
        st.markdown("---")
        st.header("Lịch sử")
        if st.session_state.history:
            for i, item in enumerate(reversed(st.session_state.history)):
                label = item.get("title", str(item)) if isinstance(item, dict) else str(item)
                col1, col2 = st.columns([5, 1])
                with col1:
                    st.button(label, key=f"history_btn_{i}", use_container_width=True,
                              on_click=handle_sidebar_click, args=(label,))
                with col2:
                    if is_logged_in:
                        if st.button("🗑️", key=f"delete_history_btn_{i}", use_container_width=True):
                            success, message = firebase.delete_from_history(st.session_state.firebase_db, st.session_state.user_info, label)
                            st.toast(message)
                            if success:
                                st.session_state.history.remove(item)
                                st.rerun()
        else:
            st.info("Chưa có lịch sử tra cứu nào.")

        # --- Bộ sưu tập ---
        if is_logged_in:
            st.markdown("---")
            st.header("Bộ sưu tập")

            # Tạo bộ sưu tập mới
            def handle_create_collection():
                coll_name = st.session_state.new_collection_input
                if coll_name:
                    success, message = firebase.create_new_collection(st.session_state.firebase_db, st.session_state.user_info, coll_name)
                    if success:
                        st.toast(message)
                        st.session_state.collections[coll_name] = []
                        st.session_state.new_collection_input = ""
                        st.rerun()
                    else:
                        st.error(message)

            st.text_input("Tên bộ sưu tập mới:", key="new_collection_input")
            st.button("Tạo mới", on_click=handle_create_collection)

            collections = st.session_state.get("collections", {})
            for name, drugs_or_placeholder in list(collections.items()):
                drugs = drugs_or_placeholder if isinstance(drugs_or_placeholder, list) else []
                expander_title = f"{name} ({len(drugs)})"
                with st.expander(expander_title):
                    if not drugs:
                        st.write("Bộ sưu tập này trống.")
                    else:
                        for idx, drug in enumerate(drugs):
                            col1, col2 = st.columns([5,1])
                            with col1:
                                st.button(drug, key=f"collection_btn_{name}_{idx}", use_container_width=True,
                                          on_click=handle_sidebar_click, args=(drug,))
                            with col2:
                                if st.button("🗑️", key=f"delete_drug_{name}_{idx}", use_container_width=True):
                                    success, message = firebase.delete_from_collection(st.session_state.firebase_db, st.session_state.user_info, name, drug)
                                    st.toast(message)
                                    if success:
                                        st.session_state.collections[name].remove(drug)
                                        st.rerun()

                    # Nút xóa bộ sưu tập
                    if st.session_state.get("confirming_delete_collection") == name:
                        if st.button(f"🔴 Xác nhận xóa '{name}'", key=f"confirm_delete_{name}"):
                            success, message = firebase.delete_collection(st.session_state.firebase_db, st.session_state.user_info, name)
                            st.toast(message)
                            if success:
                                del st.session_state.collections[name]
                            st.session_state.confirming_delete_collection = None
                            st.rerun()
                    else:
                        if st.button("🗑️", key=f"delete_collection_btn_{name}", use_container_width=True):
                            st.session_state.confirming_delete_collection = name
                            st.rerun()

        # --- Pro Access ---
        st.markdown("---")
        st.header("Pro Access")
        if st.session_state.get("pro_access"):
            st.success("Bạn đã có quyền truy cập Pro 🎉")
        else:
            if st.button("Nâng cấp Pro"):
                st.session_state.pro_access = True
                st.toast("Đã nâng cấp Pro!")
