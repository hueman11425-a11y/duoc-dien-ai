import streamlit as st
import google.generativeai as genai
import google.generativeai.types as genai_types
import google.api_core.exceptions as ga_ex
import datetime
from datetime import date, timedelta
import pandas as pd
import gspread
from gspread_dataframe import get_as_dataframe
from gspread.exceptions import SpreadsheetNotFound
from Bio import Entrez
import time
import pyrebase # Thư viện mới

# --- CẤU HÌNH FIREBASE ---
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
    auth = firebase.auth()
    st.session_state.firebase_auth = auth
except Exception as e:
    st.error("Lỗi khi khởi tạo Firebase. Vui lòng kiểm tra file secrets.toml của bạn.")
    st.stop()


# --- KIỂM TRA TRẠNG THÁI BẢO TRÌ ---
# ... (Giữ nguyên)

# --- 1. KHỞI TẠO TRẠNG THÁI PHIÊN ---
if 'history' not in st.session_state: st.session_state.history = []
if 'pro_access' not in st.session_state: st.session_state.pro_access = False
if 'user_info' not in st.session_state: st.session_state.user_info = None

# --- 2. CÁC HÀM XỬ LÝ (Giữ nguyên các hàm cũ) ---
# ... (Toàn bộ các hàm xử lý cũ của bạn từ get_access_codes_df đến run_lookup được giữ nguyên ở đây) ...

# --- 5. GIAO DIỆN VÀ LOGIC CHÍNH ---
st.set_page_config(page_title="Dược Điển AI", page_icon="💊")
st.title("Dược Điển AI 💊")
st.caption("Dự án được phát triển bởi group CÂCK và AI")

# --- HỆ THỐNG ĐĂNG NHẬP MỚI ---
auth = st.session_state.firebase_auth

# Nếu người dùng đã đăng nhập, chào mừng họ
if st.session_state.user_info:
    user_email = st.session_state.user_info['email']
    st.sidebar.success(f"Chào mừng, {user_email}")
    if st.sidebar.button("Đăng xuất"):
        st.session_state.user_info = None
        st.rerun()

    # --- LUỒNG 1: GIAO DIỆN CHÍNH KHI ĐÃ ĐĂNG NHẬP ---
    # ... (Toàn bộ giao diện chính của bạn khi đã đăng nhập được đặt ở đây)

# Nếu người dùng chưa đăng nhập, hiển thị form
else:
    choice = st.sidebar.selectbox("Đăng nhập / Đăng ký", ["Đăng nhập", "Đăng ký", "Tiếp tục với tư cách khách"])

    if choice == "Đăng nhập":
        with st.sidebar.form("login_form"):
            email = st.text_input("Email")
            password = st.text_input("Mật khẩu", type="password")
            login_button = st.form_submit_button("Đăng nhập")
            if login_button:
                try:
                    user = auth.sign_in_with_email_and_password(email, password)
                    st.session_state.user_info = user
                    st.rerun()
                except Exception as e:
                    st.sidebar.error("Email hoặc mật khẩu không chính xác.")

    elif choice == "Đăng ký":
        with st.sidebar.form("register_form"):
            email = st.text_input("Email")
            password = st.text_input("Mật khẩu", type="password")
            register_button = st.form_submit_button("Đăng ký")
            if register_button:
                try:
                    user = auth.create_user_with_email_and_password(email, password)
                    st.sidebar.success("Đăng ký thành công! Vui lòng chuyển qua tab 'Đăng nhập'.")
                except Exception as e:
                    st.sidebar.error("Email này có thể đã tồn tại hoặc không hợp lệ.")

    # --- LUỒNG 2: DÀNH CHO KHÁCH ---
    else: # choice == "Tiếp tục với tư cách khách"
        st.info("Bạn đang sử dụng với tư cách khách. Một số tính năng sẽ bị hạn chế.")
        # ... (Giao diện tra cứu cơ bản cho khách được đặt ở đây)
