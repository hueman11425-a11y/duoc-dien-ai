import streamlit as st
import pyrebase

@st.cache_resource
def initialize_firebase_app():
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
        st.error("Lỗi khi khởi tạo Firebase. Vui lòng kiểm tra file secrets.toml.")
        st.exception(e)
        return None

def display_auth_forms(auth, db):
    """
    Hiển thị các form đăng nhập/đăng ký trong sidebar.
    Trả về True nếu người dùng đã đăng nhập.
    """
    if 'user_info' not in st.session_state:
        st.session_state.user_info = None

    if st.session_state.user_info:
        user_email = st.session_state.user_info['email']
        st.success(f"Chào mừng, {user_email}")
        if st.button("Đăng xuất"):
            # Giữ lại các biến của firebase, xóa các biến người dùng
            for key in list(st.session_state.keys()):
                if key not in ['firebase_app', 'firebase_auth', 'firebase_db']:
                    del st.session_state[key]
            st.rerun()
    else:
        choice = st.selectbox("Đăng nhập / Đăng ký", ["Tiếp tục với tư cách khách", "Đăng nhập", "Đăng ký"])
        if choice == "Đăng nhập":
            with st.form("login_form"):
                email = st.text_input("Email")
                password = st.text_input("Mật khẩu", type="password")
                if st.form_submit_button("Đăng nhập"):
                    try:
                        user = auth.sign_in_with_email_and_password(email, password)
                        st.session_state.user_info = user
                        st.rerun()
                    except Exception:
                        st.error("Email hoặc mật khẩu không chính xác.")
        elif choice == "Đăng ký":
            with st.form("register_form"):
                email = st.text_input("Email")
                password = st.text_input("Mật khẩu", type="password")
                if st.form_submit_button("Đăng ký"):
                    try:
                        user = auth.create_user_with_email_and_password(email, password)
                        st.sidebar.success("Đăng ký thành công! Đang khởi tạo dữ liệu...")
                        user_id = user['localId']
                        token = user['idToken']
                        default_data = {"history": [], "collections": {}, "is_pro": False, "usage_counters": {"prescription_analysis": 0}}
                        db.child("user_data").child(user_id).set(default_data, token=token)
                        st.info("Vui lòng chuyển qua tab 'Đăng nhập'.")
                    except Exception:
                        st.sidebar.error("Email này có thể đã tồn tại hoặc không hợp lệ.")

    return st.session_state.user_info is not None
