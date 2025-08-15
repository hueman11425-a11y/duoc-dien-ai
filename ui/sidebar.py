import streamlit as st

def render_sidebar(verify_code):
    st.sidebar.header("Lịch sử tra cứu")
    if not st.session_state.history:
        st.sidebar.info("Chưa có thuốc nào được tra cứu.")
    else:
        for drug in st.session_state.history:
            if st.sidebar.button(drug, key=f"history_{drug}", use_container_width=True):
                st.session_state.selected_drug = drug
                st.session_state.trigger_lookup = True

    st.sidebar.markdown("---")

    # Giữ nguyên đúng bố cục bản gốc: container có viền
    with st.sidebar.container(border=True):
        st.write("**Bạn có ý tưởng để cải thiện ứng dụng?**")
        st.link_button(
            "Gửi phản hồi!",
            url="https://forms.gle/M44GDS4hJ7LpY7b98",
            help="Mở form góp ý trong một tab mới"
        )

    st.sidebar.markdown("---")

    st.sidebar.header("Truy cập Pro")
    if st.session_state.get("pro_access"):
        st.sidebar.success("Bạn đã có quyền truy cập Pro.")
    else:
        code = st.sidebar.text_input("Nhập mã Pro:", type="password")
        if st.sidebar.button("Xác thực"):
            is_valid, message = verify_code(code)
            if is_valid:
                st.sidebar.success(message)
                st.rerun()
            else:
                st.sidebar.error(message)
