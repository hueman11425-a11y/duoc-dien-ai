import streamlit as st

def init_session_state():
    defaults = {
        "current_page": "Tra cứu Dược điển",
        "user_info": None,
        "user_data_loaded": False,
        "history": [],
        "collections": {},
        "pro_access": False,
        "confirming_delete_collection": None,
        "guest_cache": {},
        "query_result": None,
        "analysis_result": None,
        "lookup_input_trigger": None,
    }
    for k, v in defaults.items():
        if k not in st.session_state:
            st.session_state[k] = v
