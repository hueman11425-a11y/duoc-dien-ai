import streamlit as st
import gspread
import pandas as pd
from gspread_dataframe import get_as_dataframe

@st.cache_data(ttl=600)
def get_access_codes_df():
    try:
        credentials_dict = dict(st.secrets.connections.gsheets.credentials)
        gspread_client = gspread.service_account_from_dict(credentials_dict)
        
        spreadsheet_name = st.secrets.connections.gsheets.spreadsheet
        spreadsheet = gspread_client.open(spreadsheet_name)
        worksheet = spreadsheet.sheet1
        return get_as_dataframe(worksheet)
    except Exception:
        st.error("Lỗi kết nối tới Google Sheets.")
        return pd.DataFrame()

def verify_code(db, user_info, user_code):
    from utils.constants import COLLECTION_LIMIT
    if not user_code:
        return False, "Vui lòng nhập mã truy cập."

    user_id = user_info['localId']
    token = user_info['idToken']
    codes_df = get_access_codes_df()
    if codes_df.empty:
        return False, "Không thể tải dữ liệu mã truy cập."
    
    codes_df.dropna(subset=['code'], inplace=True)
    codes_df['code'] = codes_df['code'].astype(str)

    if user_code not in codes_df['code'].values:
        return False, "Mã không hợp lệ hoặc không tìm thấy."
    
    claimed_codes = db.child("claimed_pro_codes").get(token=token).val() or {}
    if user_code in claimed_codes:
        if claimed_codes[user_code] != user_id:
            return False, "Mã này đã được sử dụng bởi một tài khoản khác."
        else:
            st.session_state.pro_access = True
            return True, "Bạn đã kích hoạt mã này trước đó. Chào mừng trở lại!"
    try:
        db.child("claimed_pro_codes").child(user_code).set(user_id, token=token)
        db.child("user_data").child(user_id).child("is_pro").set(True, token=token)
        st.session_state.pro_access = True
        return True, "Kích hoạt PRO thành công! Cảm ơn bạn đã ủng hộ dự án."
    except Exception as e:
        return False, f"Đã xảy ra lỗi trong quá trình kích hoạt: {e}"
