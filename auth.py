import streamlit as st
import pyrebase
import logging

# Cấu hình logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

@st.cache_resource
def initialize_firebase_app():
    """
    Khởi tạo kết nối tới Firebase và trả về đối tượng app chính.
    Sử dụng cache để tránh khởi tạo lại nhiều lần khi rerun.
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
        logging.info("✅ Firebase initialized successfully")
        return firebase
    except Exception as e:
        logging.error("❌ Firebase initialization failed", exc_info=True)
        st.error("Lỗi khi khởi tạo Firebase. Vui lòng kiểm tra `secrets.toml` và cấu hình Firebase.")
        return None


def reset_user_session():
    """Xóa toàn bộ thông tin liên quan đến người dùng trong session_state."""
    keys_to_reset = [
        "user_info", "history", "pro_access",
        "user_data_loaded", "collections", "last_drug_searched"
    ]
    for key in keys_to_reset:
        st.session_state[key] = None if key == "user_info" else []
    st.session_state.collections = {}
    logging.info("🔄 User session reset")


def display_auth_forms(auth):
    """
    Hiển thị form xác thực người dùng (login, register, guest).
    Trả về True nếu có người dùng đăng nhập, False nếu không.
    """
    if "user_info" not in st.session_state:
        st.session_state.user_info = None

    with st.sidebar:
        # Trường hợp đã đăng nhập
        if st.session_state.user_info:
            user_email = st.session_state.user_info.get("email", "Người dùng")
            st.success(f"👋 Chào mừng, {user_email}")
            if st.button("Đăng xuất"):
                reset_user_session()
                st.rerun()
        
        # Trường hợp chưa đăng nhập
        else:
            choice = st.selectbox(
                "Đăng nhập / Đăng ký",
                ["Tiếp tục với tư cách khách", "Đăng nhập", "Đăng ký"]
            )

            if choice == "Đăng nhập":
                with st.form("login_form"):
                    email = st.text_input("Email")
                    password = st.text_input("Mật khẩu", type="password")
                    if st.form_submit_button("Đăng nhập"):
                        try:
                            user = auth.sign_in_with_email_and_password(email, password)
                            st.session_state.user_info = user
                            logging.info(f"✅ User logged in: {email}")
                            st.rerun()
                        except Exception as e:
                            logging.warning(f"❌ Login failed for {email}", exc_info=True)
                            st.error("Email hoặc mật khẩu không chính xác.")

            elif choice == "Đăng ký":
                with st.form("register_form"):
                    email = st.text_input("Email")
                    password = st.text_input("Mật khẩu", type="password")
                    if st.form_submit_button("Đăng ký"):
                        try:
                            auth.create_user_with_email_and_password(email, password)
                            st.sidebar.success("🎉 Đăng ký thành công! Hãy đăng nhập ngay.")
                            logging.info(f"✅ User registered: {email}")
                        except Exception as e:
                            error_str = str(e)
                            logging.error(f"❌ Registration error: {error_str}", exc_info=True)

                            if "OPERATION_NOT_ALLOWED" in error_str:
                                st.sidebar.error("❌ Firebase Authentication chưa bật Email/Password.\n👉 Vào Firebase Console → Authentication → Sign-in method để bật.")
                            elif "EMAIL_EXISTS" in error_str:
                                st.sidebar.error("❌ Email đã tồn tại. Hãy dùng email khác.")
                            elif "WEAK_PASSWORD" in error_str:
                                st.sidebar.error("❌ Mật khẩu quá yếu (ít nhất 6 ký tự).")
                            elif "INVALID_EMAIL" in error_str:
                                st.sidebar.error("❌ Email không hợp lệ. Vui lòng nhập đúng định dạng.")
                            else:
                                st.sidebar.error("❌ Lỗi không xác định khi đăng ký. Xem log để biết thêm chi tiết.")

    return st.session_state.user_info is not None
