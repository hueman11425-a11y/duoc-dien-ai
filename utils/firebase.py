import streamlit as st
from utils.constants import HISTORY_LIMIT, COLLECTION_LIMIT, DRUGS_PER_COLLECTION_LIMIT
from html import unescape

def load_user_data(db, user_info):
    try:
        user_id = user_info['localId']
        token = user_info['idToken']
        data = db.child("user_data").child(user_id).get(token=token).val()
        if not data:
            return [], {}  # chỉ 2 giá trị
        history = data.get("history", [])
        collections = data.get("collections", {})
        is_pro = data.get("is_pro", False)
        if is_pro:
            st.session_state.pro_access = True
        return history, collections
    except Exception:
        return [], {}
    

def load_user_result(db, user_info, drug_name):
    try:
        user_id = user_info['localId']
        token = user_info['idToken']
        result = db.child("user_data").child(user_id).child("results_cache").child(drug_name).get(token=token).val()
        return result
    except Exception:
        return None

def save_new_result(db, user_info, drug_name, result_text):
    try:
        user_id = user_info['localId']
        token = user_info['idToken']
        is_pro = db.child("user_data").child(user_id).child("is_pro").get(token=token).val() or False
        db.child("user_data").child(user_id).child("results_cache").child(drug_name).set(result_text, token=token)

        history = db.child("user_data").child(user_id).child("history").get(token=token).val() or []
        if drug_name in history:
            history.remove(drug_name)
        history.insert(0, drug_name)

        if not is_pro and len(history) > HISTORY_LIMIT:
            drug_to_delete = history.pop()
            db.child("user_data").child(user_id).child("results_cache").child(drug_to_delete).remove(token=token)

        db.child("user_data").child(user_id).child("history").set(history, token=token)
        return history
    except Exception as e:
        st.error(f"Lỗi khi lưu kết quả tra cứu: {e}")
        return None

def create_new_collection(db, user_info, collection_name):
    if not collection_name or collection_name.isspace():
        return False, "Tên bộ sưu tập không được để trống."
    try:
        user_id = user_info['localId']
        token = user_info['idToken']
        is_pro = db.child("user_data").child(user_id).child("is_pro").get(token=token).val() or False
        collections = db.child("user_data").child(user_id).child("collections").get(token=token).val() or {}

        if not is_pro and len(collections) >= COLLECTION_LIMIT and collection_name not in collections:
            return False, f"Đã đạt giới hạn {COLLECTION_LIMIT} bộ sưu tập."
        if collection_name in collections:
            return False, f"Bộ sưu tập '{collection_name}' đã tồn tại."

        collections[collection_name] = True
        db.child("user_data").child(user_id).child("collections").set(collections, token=token)
        return True, f"Đã tạo thành công bộ sưu tập '{collection_name}'."
    except Exception as e:
        return False, f"Đã xảy ra lỗi không xác định: {e}"

def add_drug_to_collection(db, user_info, collection_name, drug_name):
    try:
        user_id = user_info['localId']
        token = user_info['idToken']
        is_pro = db.child("user_data").child(user_id).child("is_pro").get(token=token).val() or False
        collections = db.child("user_data").child(user_id).child("collections").get(token=token).val() or {}

        if collection_name not in collections:
            return f"Lỗi: Không tìm thấy bộ sưu tập '{collection_name}'."

        drug_list_or_placeholder = collections.get(collection_name)
        if drug_list_or_placeholder is True:
            drug_list = []
        elif isinstance(drug_list_or_placeholder, list):
            drug_list = drug_list_or_placeholder
        else:
            return "Lỗi: Cấu trúc dữ liệu của bộ sưu tập không hợp lệ."

        if drug_name in drug_list:
            return f"'{drug_name}' đã có trong bộ sưu tập này."
        if not is_pro and len(drug_list) >= DRUGS_PER_COLLECTION_LIMIT:
            return f"Bộ sưu tập '{collection_name}' đã đầy (tối đa {DRUGS_PER_COLLECTION_LIMIT} thuốc)."


        drug_name_clean = unescape(drug_name).strip()
        if not drug_name_clean:
            return "Tên thuốc không hợp lệ."
        if drug_name_clean in drug_list:
            return f"'{drug_name_clean}' đã có trong bộ sưu tập này."
        drug_list.append(drug_name_clean)
        collections[collection_name] = drug_list
        db.child("user_data").child(user_id).child("collections").set(collections, token=token)
        return f"Đã thêm '{drug_name}' vào '{collection_name}'."
   
    except Exception as e:
        return f"Đã xảy ra lỗi không xác định: {e}"

def delete_from_history(db, user_info, drug_name):
    try:
        user_id = user_info['localId']
        token = user_info['idToken']
        history = db.child("user_data").child(user_id).child("history").get(token=token).val() or []
        if drug_name in history:
            history.remove(drug_name)
            db.child("user_data").child(user_id).child("history").set(history, token=token)
        return True, f"Đã xóa '{drug_name}' khỏi lịch sử."
    except Exception as e:
        return False, f"Lỗi khi xóa khỏi lịch sử: {e}"

def delete_from_collection(db, user_info, collection_name, drug_name):
    try:
        user_id = user_info['localId']
        token = user_info['idToken']
        collections = db.child("user_data").child(user_id).child("collections").get(token=token).val() or {}
        if collection_name in collections and isinstance(collections[collection_name], list):
            drug_list = collections[collection_name]
            if drug_name in drug_list:
                drug_list.remove(drug_name)
                if not drug_list:
                    collections[collection_name] = True
                else:
                    collections[collection_name] = drug_list
                db.child("user_data").child(user_id).child("collections").set(collections, token=token)
                return True, f"Đã xóa '{drug_name}' khỏi '{collection_name}'."
        return False, "Không tìm thấy thuốc hoặc bộ sưu tập."
    except Exception as e:
        return False, f"Lỗi khi xóa khỏi bộ sưu tập: {e}"

def delete_collection(db, user_info, collection_name):
    try:
        user_id = user_info['localId']
        token = user_info['idToken']
        db.child("user_data").child(user_id).child("collections").child(collection_name).remove(token=token)
        return True, f"Đã xóa bộ sưu tập '{collection_name}'."
    except Exception as e:
        return False, f"Lỗi khi xóa bộ sưu tập: {e}"
