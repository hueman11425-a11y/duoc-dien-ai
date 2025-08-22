import streamlit as st
import google.generativeai as genai

@st.cache_resource
def get_lookup_model():
    model_name = st.secrets.get("models", {}).get("lookup", "gemini-1.5-flash-latest")
    return genai.GenerativeModel(model_name)

@st.cache_resource
def get_pro_model():
    model_name = st.secrets.get("models", {}).get("pro", "gemini-pro")
    return genai.GenerativeModel(model_name)

@st.cache_resource
def get_prescription_model():
    # model_name = st.secrets.get("models", {}).get("prescription", "gemini-1.5-pro-latest")
    model_name = st.secrets.get("models", {}).get("prescription", "gemini-1.5-flash-latest")
    return genai.GenerativeModel(model_name)
