import streamlit as st
from datetime import date
from utils.models import get_prescription_model
from utils.prompts import PROMPT_PRESCRIPTION
from utils.constants import PRESCRIPTION_LIMIT_PER_DAY

def get_prescription_analysis(db, user_info, patient_context, prescription_text):
    try:
        user_id = user_info['localId']
        token = user_info['idToken']
        is_pro = st.session_state.get("pro_access", False)
        
        if not is_pro:
            today_str = date.today().isoformat()
            usage_ref = db.child("user_data").child(user_id).child("usage_counters").child("prescription_analysis")
            usage_data = usage_ref.get(token=token).val() or {"count": 0, "last_updated": ""}
            
            if usage_data.get("last_updated") != today_str:
                usage_data = {"count": 0, "last_updated": today_str}
            
            if usage_data.get("count", 0) >= PRESCRIPTION_LIMIT_PER_DAY:
                return f"❌ Bạn đã hết {PRESCRIPTION_LIMIT_PER_DAY} lượt phân tích miễn phí trong ngày hôm nay."

        model = get_prescription_model()
        prompt = PROMPT_PRESCRIPTION.format(patient_context=patient_context, prescription_text=prescription_text)
        response = model.generate_content(prompt)
        
        if not is_pro:
            usage_data["count"] += 1
            usage_ref.set(usage_data, token=token)
            
        return response.text
    except Exception as e:
        return f" Lỗi trong quá trình phân tích: {e}"
