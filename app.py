import streamlit as st
import auth
import utils

# --- KHỞI TẠO LẦN ĐẦU ---
st.set_page_config(page_title="Dược Điển AI", page_icon="💊", layout="wide")

# Khởi tạo các dịch vụ
if "firebase_app" not in st.session_state:
    st.session_state.firebase_app = auth.initialize_firebase_app()
    if st.session_state.firebase_app:
        st.session_state.firebase_auth = st.session_state.firebase_app.auth()
        st.session_state.firebase_db = st.session_state.firebase_app.database()
    else:
        st.stop()

# Khởi tạo các biến trạng thái
if 'current_page' not in st.session_state: st.session_state.current_page = "Tra cứu Dược điển"
if 'user_info' not in st.session_state: st.session_state.user_info = None
if 'user_data_loaded' not in st.session_state: st.session_state.user_data_loaded = False
if 'history' not in st.session_state: st.session_state.history = []
if 'collections' not in st.session_state: st.session_state.collections = {}
if 'pro_access' not in st.session_state: st.session_state.pro_access = False
if 'confirming_delete_collection' not in st.session_state: st.session_state.confirming_delete_collection = None
if 'guest_cache' not in st.session_state: st.session_state.guest_cache = {}
if 'query_result' not in st.session_state: st.session_state.query_result = None
if 'analysis_result' not in st.session_state: st.session_state.analysis_result = None
if 'lookup_input_trigger' not in st.session_state: st.session_state.lookup_input_trigger = None


# --- HÀM CALLBACK CHO LỊCH SỬ VÀ BỘ SƯU TẬP ---
def handle_sidebar_click(drug_name):
    """
    Callback an toàn để xử lý khi một thuốc được chọn từ sidebar (Lịch sử hoặc Bộ sưu tập).
    Nó sẽ tải kết quả đã lưu thay vì thực hiện một tra cứu mới.
    """
    st.session_state.current_page = "Tra cứu Dược điển"
    st.session_state.lookup_input_field = drug_name # Cập nhật ô text input
    
    is_logged_in = st.session_state.user_info is not None
    if is_logged_in:
        result = utils.load_user_result(st.session_state.firebase_db, st.session_state.user_info, drug_name)
    else:
        result = st.session_state.guest_cache.get(drug_name)
    
    st.session_state.query_result = result if result else "Không tìm thấy kết quả đã lưu."

# --- CÁC HÀM ĐỂ VẼ NỘI DUNG TỪNG TRANG ---
def render_lookup_page():
    st.header("Tra cứu Dược điển 💊")
    
    is_logged_in = st.session_state.user_info is not None
    firebase_db = st.session_state.firebase_db
    user_info = st.session_state.user_info
    
    def run_lookup(drug_name):
        st.session_state.query_result = None
        is_pro = st.session_state.get("pro_access", False)
        with st.spinner(f"Đang tra cứu '{drug_name}'..."):
            api_result, identified_name = utils.get_drug_info_from_api(drug_name, is_pro)
            st.session_state.query_result = api_result
            if identified_name and "không được nhận dạng" not in api_result:
                if is_logged_in:
                    updated_history = utils.save_new_result(firebase_db, user_info, identified_name, api_result)
                    if updated_history is not None: st.session_state.history = updated_history
                else:
                    st.session_state.guest_cache[identified_name] = api_result
                    if identified_name not in st.session_state.history:
                        st.session_state.history.insert(0, identified_name)
                        if len(st.session_state.history) > 10:
                            drug_to_remove = st.session_state.history.pop()
                            if drug_to_remove in st.session_state.guest_cache:
                                del st.session_state.guest_cache[drug_to_remove]

    search_query = st.text_input("Nhập tên thuốc (biệt dược hoặc hoạt chất):", key="lookup_input_field")
    
    if st.button("Tra cứu"):
        if search_query:
            run_lookup(search_query)
        else:
            st.warning("Vui lòng nhập tên thuốc trước khi tra cứu.")

    if st.session_state.query_result:
        result_to_display = st.session_state.query_result
        if result_to_display.startswith("❌ Lỗi:"):
            st.error(result_to_display)
        else:
            st.markdown(result_to_display)
            if is_logged_in:
                identified_name = result_to_display.split("**")[1] if "**" in result_to_display else None
                if identified_name:
                    st.markdown("---")
                    st.subheader(f"Lưu '{identified_name}' vào bộ sưu tập")
                    collections = st.session_state.get("collections", {})
                    if not collections:
                        st.info("Bạn chưa có bộ sưu tập nào.")
                    else:
                        col1, col2 = st.columns([2,1])
                        with col1:
                            sel_coll = st.selectbox("Chọn bộ sưu tập:", options=list(collections.keys()))
                        with col2:
                            st.write(""); st.write("")
                            if st.button("Thêm thuốc", use_container_width=True):
                                msg = utils.add_drug_to_collection(firebase_db, user_info, sel_coll, identified_name)
                                st.toast(msg)
                                _, collections_new = utils.load_user_data(firebase_db, user_info)
                                st.session_state.collections = collections_new
                                st.rerun()

def render_prescription_analysis_page():
    st.header("Phân tích Đơn thuốc 🩺")
    is_logged_in = st.session_state.user_info is not None
    
    if not is_logged_in:
        st.warning("Vui lòng đăng nhập để sử dụng tính năng này.")
        return
    
    st.subheader("Phần 1: Bối cảnh Bệnh nhân (Không bắt buộc)")
    col1, col2 = st.columns(2)
    with col1:
        conditions = st.text_area("Tình trạng bệnh lý nền:", placeholder="Ví dụ: Suy thận mạn, tăng huyết áp...", height=150)
    with col2:
        allergies = st.text_area("Dị ứng thuốc đã biết:", placeholder="Ví dụ: Dị ứng với penicillin...", height=150)
        
    st.subheader("Phần 2: Thông tin Đơn thuốc (Bắt buộc)")
    prescription_text = st.text_area("Dán nội dung đơn thuốc vào đây:", placeholder="1. Tên thuốc A 50mg - 1 viên x 2...", height=250)
    
    if st.button("Phân tích Đơn thuốc", type="primary"):
        if not prescription_text.strip():
            st.error("Vui lòng nhập thông tin đơn thuốc.")
        else:
            with st.spinner("AI đang phân tích, vui lòng chờ trong giây lát..."):
                patient_context = f"- Bệnh lý nền: {conditions if conditions.strip() else 'Không có'}\n- Dị ứng: {allergies if allergies.strip() else 'Không có'}"
                result = utils.get_prescription_analysis(st.session_state.firebase_db, st.session_state.user_info, patient_context, prescription_text)
                st.session_state.analysis_result = result
    
    if st.session_state.analysis_result:
        st.markdown("---")
        st.subheader("Kết quả Phân tích")
        if st.session_state.analysis_result.startswith("❌"):
             st.error(st.session_state.analysis_result)
        else:
            st.markdown(st.session_state.analysis_result)

# --- SIDEBAR VÀ ĐIỀU HƯỚNG ---
with st.sidebar:
    is_logged_in = auth.display_auth_forms(st.session_state.firebase_auth, st.session_state.firebase_db)
    
    if is_logged_in and not st.session_state.user_data_loaded:
        with st.spinner("Đang tải dữ liệu..."):
            history, collections = utils.load_user_data(st.session_state.firebase_db, st.session_state.user_info)
            st.session_state.history = history
            st.session_state.collections = collections
        st.session_state.user_data_loaded = True
        
    st.markdown("---")
    
    st.header("Tính năng")
    page_options = ["Tra cứu Dược điển", "Phân tích Đơn thuốc"]
    
    def on_page_change():
        st.session_state.current_page = st.session_state.feature_radio
        st.session_state.query_result = None
        st.session_state.analysis_result = None

    try:
        current_page_index = page_options.index(st.session_state.current_page)
    except ValueError:
        current_page_index = 0

    st.radio(
        "Chọn tính năng:", page_options, 
        index=current_page_index,
        on_change=on_page_change, 
        key='feature_radio',
        label_visibility="collapsed"
    )
    
    st.markdown("---")
    
    st.header("Lịch sử tra cứu")
    history_container = st.container(height=300)
    with history_container:
        if not st.session_state.history:
            st.info("Chưa có thuốc nào được tra cứu.")
        else:
            for drug in st.session_state.history[:]:
                col1, col2 = st.columns([5, 1])
                with col1:
                    st.button(drug, key=f"history_{drug}", use_container_width=True, on_click=handle_sidebar_click, args=(drug,))
                with col2:
                    if is_logged_in:
                        if st.button("🗑️", key=f"delete_history_{drug}", use_container_width=True):
                            success, message = utils.delete_from_history(st.session_state.firebase_db, st.session_state.user_info, drug)
                            st.toast(message)
                            if success:
                                st.session_state.history.remove(drug)
                                st.rerun()
    st.markdown("---")
    if is_logged_in:
        st.header("Bộ sưu tập")
        def handle_create_collection():
            coll_name = st.session_state.new_collection_input
            if coll_name:
                success, message = utils.create_new_collection(st.session_state.firebase_db, st.session_state.user_info, coll_name)
                if success:
                    st.toast(message); st.session_state.collections[coll_name] = True; st.session_state.new_collection_input = ""
                else: st.error(message)
        st.text_input("Tên bộ sưu tập mới:", key="new_collection_input")
        st.button("Tạo mới", on_click=handle_create_collection)
        collections = st.session_state.get("collections", {})
        is_pro = st.session_state.get("pro_access", False)
        for name, drugs_or_placeholder in list(collections.items()):
            coll_col1, coll_col2 = st.columns([5,1])
            with coll_col1:
                if st.session_state.confirming_delete_collection == name:
                    if st.button(f"🔴 Xác nhận xóa '{name}'", key=f"confirm_delete_{name}", use_container_width=True):
                        success, message = utils.delete_collection(st.session_state.firebase_db, st.session_state.user_info, name)
                        st.toast(message); st.session_state.confirming_delete_collection = None
                        if success: del st.session_state.collections[name]
                        st.rerun()
                else:
                    drugs = drugs_or_placeholder if isinstance(drugs_or_placeholder, list) else []
                    expander_title = f"{name} ({len(drugs)}{'' if is_pro else f'/{utils.DRUGS_PER_COLLECTION_LIMIT}'} thuốc)"
                    with st.expander(expander_title):
                        if not drugs: st.write("Bộ sưu tập này trống.")
                        else:
                            for drug in drugs:
                                d_col1, d_col2 = st.columns([5, 1])
                                with d_col1:
                                    st.button(drug, key=f"collection_{name}_{drug}", use_container_width=True, on_click=handle_sidebar_click, args=(drug,))
                                with d_col2:
                                    if st.button("🗑️", key=f"delete_collection_drug_{name}_{drug}", use_container_width=True):
                                        success, message = utils.delete_from_collection(st.session_state.firebase_db, st.session_state.user_info, name, drug)
                                        st.toast(message)
                                        if success:
                                            _, collections_new = utils.load_user_data(st.session_state.firebase_db, st.session_state.user_info)
                                            st.session_state.collections = collections_new
                                            st.rerun()
            with coll_col2:
                if st.session_state.confirming_delete_collection != name:
                    if st.button("🗑️", key=f"delete_collection_{name}", use_container_width=True):
                        st.session_state.confirming_delete_collection = name
                        st.rerun()
        if is_pro: st.markdown(f"Đã tạo {len(collections)} bộ sưu tập (PRO).")
        else: st.markdown(f"Đã tạo {len(collections)}/{utils.COLLECTION_LIMIT} bộ sưu tập.")
        st.markdown("---")
    with st.container(border=True):
        st.write("**Bạn có ý tưởng để cải thiện ứng dụng?**")
        st.link_button("Gửi phản hồi ngay!", url="https://forms.gle/M44GDS4hJ7LpY7b98")
    if is_logged_in:
        st.header("Truy cập Pro")
        if st.session_state.get("pro_access"): st.success("Bạn đã có quyền truy cập Pro.")
        else:
            pro_code_input = st.text_input("Nhập mã truy cập Pro:", type="password")
            if st.button("Xác thực"):
                user_info = st.session_state.get("user_info")
                if user_info:
                    is_valid, message = utils.verify_code(st.session_state.firebase_db, user_info, pro_code_input)
                    if is_valid: st.success(message); st.rerun()
                    else: st.error(message)
                else: st.warning("Vui lòng đăng nhập để kích hoạt mã.")

# --- HIỂN THỊ NỘI DUNG TRANG CHÍNH DỰA TRÊN LỰA CHỌN ---
st.title("Dược Điển AI")
if st.session_state.current_page == "Tra cứu Dược điển":
    render_lookup_page()
elif st.session_state.current_page == "Phân tích Đơn thuốc":
    render_prescription_analysis_page()
