# Xác thực mã Pro
import pandas as pd
import gspread
from gspread_dataframe import get_as_dataframe
from datetime import date
import streamlit as st

@st.cache_data(ttl=600)
def get_access_codes_df() -> pd.DataFrame:
    try:
        credentials = st.secrets.connections.gsheets.credentials
        client = gspread.service_account_from_dict(credentials)
        sheet = client.open(st.secrets.connections.gsheets.spreadsheet).sheet1
        return get_as_dataframe(sheet)
    except Exception:
        return pd.DataFrame()

def verify_code(user_code: str):
    if not user_code:
        return False, "Vui lòng nhập mã truy cập."
    
    codes_df = get_access_codes_df()
    if codes_df.empty:
        return False, "Không thể tải dữ liệu mã truy cập."
    
    codes_df.dropna(subset=['code'], inplace=True)
    codes_df['code'] = codes_df['code'].astype(str)
    match = codes_df[codes_df['code'] == user_code]
    if match.empty:
        return False, "Mã không hợp lệ."
    
    info = match.iloc[0]
    if info['type'] == 'permanent':
        st.session_state.pro_access = True
        return True, f"Chào mừng {info.get('owner', 'Pro User')}!"
    
    if info['type'] == 'temporary':
        try:
            days_passed = (date.today() - pd.to_datetime(info['created_at']).date()).days
            if 0 <= days_passed <= 7:
                st.session_state.pro_access = True
                return True, f"Mã còn hiệu lực {7 - days_passed} ngày."
            return False, "Mã tạm thời đã hết hạn."
        except Exception:
            return False, "Lỗi định dạng ngày tháng."
    
    return False, "Loại mã không xác định."
