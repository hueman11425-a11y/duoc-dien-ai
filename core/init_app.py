import streamlit as st
import auth

def init_app():
    if "firebase_app" not in st.session_state:
        st.session_state.firebase_app = auth.initialize_firebase_app()
        if st.session_state.firebase_app:
            st.session_state.firebase_auth = st.session_state.firebase_app.auth()
            st.session_state.firebase_db = st.session_state.firebase_app.database()
        else:
            st.stop()
