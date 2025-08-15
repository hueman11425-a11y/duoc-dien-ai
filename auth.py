import streamlit as st
import pyrebase

@st.cache_resource
def initialize_firebase_app():
    """
    Khởi tạo kết nối tới Firebase và trả về đối tượng app chính.
    Việc sử dụng @st.cache_resource giúp giữ kết nối ổn định qua các lần rerun.
    """
    try:
        firebase_config = {
            "apiKey": st.secrets.firebase.apiKey,
            "authDomain": st.secrets.firebase.authDomain,
            "projectId": st.secrets.firebase.projectId,
            "storageBucket": st.secrets.firebase.storageBucket,
            "messagingSenderId": st.secrets.firebase.messagingSenderId,
            "appId": st.secrets.firebase.appId,
            "databaseURL": st.secrets.firebase.databaseURL
        }
        firebase = pyrebase.initialize_app(firebase_config)
        return firebase
    except Exception as e:
        st.error("Lỗi khi khởi tạo Firebase. Vui lòng kiểm tra file secrets.toml của bạn.")
        st.exception(e)
        return None

def display_auth_forms(auth):
    """
    Hiển thị các form đăng nhập, đăng ký và xử lý logic.
    Trả về thông tin người dùng nếu đăng nhập thành công, ngược lại trả về None.
    """
    if 'user_info' not in st.session_state:
        st.session_state.user_info = None

    with st.sidebar:
        if st.session_state.user_info:
            user_email = st.session_state.user_info['email']
            st.success(f"Chào mừng, {user_email}")
            if st.button("Đăng xuất"):
                # Reset tất cả các session state liên quan đến người dùng
                st.session_state.user_info = None
                st.session_state.history = []
                st.session_state.pro_access = False
                st.session_state.user_data_loaded = False
                st.session_state.collections = {}
                st.session_state.last_drug_searched = None
                st.rerun()
        else:
            choice = st.selectbox("Đăng nhập / Đăng ký", ["Tiếp tục với tư cách khách", "Đăng nhập", "Đăng ký"])

            if choice == "Đăng nhập":
                with st.form("login_form"):
                    email = st.text_input("Email")
                    password = st.text_input("Mật khẩu", type="password")
                    login_button = st.form_submit_button("Đăng nhập")
                    if login_button:
                        try:
                            user = auth.sign_in_with_email_and_password(email, password)
                            st.session_state.user_info = user
                            st.rerun()
                        except Exception as e:
                            st.error("Email hoặc mật khẩu không chính xác.")
            
            elif choice == "Đăng ký":
                with st.form("register_form"):
                    email = st.text_input("Email")
                    password = st.text_input("Mật khẩu", type="password")
                    register_button = st.form_submit_button("Đăng ký")
                    if register_button:
                        try:
                            user = auth.create_user_with_email_and_password(email, password)
                            st.sidebar.success("Đăng ký thành công! Vui lòng chuyển qua tab 'Đăng nhập'.")
                        except Exception as e:
                            st.sidebar.error("Email này có thể đã tồn tại hoặc không hợp lệ.")

    return st.session_state.user_info is not None
