import streamlit as st
import google.generativeai as genai

# Import các file chúng ta đã tách
import auth
import utils

# --- CSS TÙY CHỈNH CHO THANH CUỘN ---
st.markdown("""
    <style>
    .scrollable-container {
        border: 1px solid #262730;
        border-radius: 0.5rem;
        padding: 10px;
        max-height: 300px; /* Chiều cao tối đa của khung */
        overflow-y: auto; /* Thêm thanh cuộn khi cần */
    }
    .scrollable-container .stButton>button {
        width: 100%;
        text-align: left;
    }
    </style>
""", unsafe_allow_html=True)


# --- KHỞI TẠO CÁC DỊCH VỤ ---
firebase_app = auth.initialize_firebase_app()
if not firebase_app:
    st.stop()

firebase_auth = firebase_app.auth()
firebase_db = firebase_app.database()

try:
    genai.configure(api_key=st.secrets["GOOGLE_API_KEY"])
except (FileNotFoundError, KeyError):
    st.error("LỖI: Vui lòng cấu hình GOOGLE_API_KEY trong secrets.toml.")
    st.stop()

# --- KHỞI TẠO TRẠNG THÁI PHIÊN ---
if 'user_data_loaded' not in st.session_state: st.session_state.user_data_loaded = False
if "query_result" not in st.session_state: st.session_state.query_result = None
# Khởi tạo history cho cả khách và người dùng đăng nhập
if 'history' not in st.session_state: st.session_state.history = []


# --- HÀM LOGIC TRUNG TÂM ---
def run_lookup(drug_name):
    st.session_state.query_result = None # Xóa kết quả cũ trước khi tra cứu
    user_info = st.session_state.get("user_info")
    is_pro = st.session_state.get("pro_access", False)

    with st.spinner(f"Đang tra cứu '{drug_name}'..."):
        # Bước 1: Kiểm tra trong kho lưu trữ vĩnh viễn (chỉ cho người dùng đăng nhập)
        if user_info:
            cached_result = utils.load_user_result(firebase_db, user_info, drug_name)
            if cached_result:
                st.session_state.query_result = cached_result
                # Cập nhật lại history list phòng trường hợp nó chưa đồng bộ
                if drug_name not in st.session_state.history:
                     st.session_state.history.insert(0, drug_name)
                return

        # Bước 2: Nếu không có, gọi API
        api_result, identified_name = utils.get_drug_info_from_api(drug_name, is_pro)
        
        # Bước 3: Xử lý kết quả
        st.session_state.query_result = api_result
        if user_info:
            # Nếu là người dùng đăng nhập, lưu kết quả mới vào kho vĩnh viễn
            if identified_name:
                updated_history = utils.save_new_result(firebase_db, user_info, identified_name, api_result)
                st.session_state.history = updated_history
        else:
            # NẾU LÀ KHÁCH, LƯU VÀO LỊCH SỬ TẠM THỜI
            if identified_name:
                if identified_name not in st.session_state.history:
                    st.session_state.history.insert(0, identified_name)
                    if len(st.session_state.history) > 10: # Giới hạn 10 mục cho khách
                        st.session_state.history.pop()


# --- BẮT ĐẦU GIAO DIỆN ---
st.set_page_config(page_title="Dược Điển AI", page_icon="💊")

# Hiển thị form đăng nhập và lấy trạng thái
is_logged_in = auth.display_auth_forms(firebase_auth)

# --- TẢI DỮ LIỆU NGƯỜI DÙNG KHI ĐĂNG NHẬP ---
if is_logged_in and not st.session_state.user_data_loaded:
    user_info = st.session_state.user_info
    with st.spinner("Đang tải dữ liệu của bạn..."):
        history, collections = utils.load_user_data(firebase_db, user_info)
        st.session_state.history = history
        st.session_state.collections = collections
    st.session_state.user_data_loaded = True

# --- GIAO DIỆN CHÍNH ---
st.title("Dược Điển AI 💊")
st.caption("Dự án được phát triển bởi group CÂCK và AI")

st.text_input("Nhập tên thuốc (biệt dược hoặc hoạt chất):", key="main_input")
if st.button("Tra cứu"):
    if st.session_state.main_input:
        run_lookup(st.session_state.main_input)
    else:
        st.warning("Vui lòng nhập tên thuốc trước khi tra cứu.")

# --- HIỂN THỊ KẾT QUẢ TRA CỨU ---
if st.session_state.query_result:
    result_to_display = st.session_state.query_result
    if result_to_display.startswith("❌ Lỗi:"):
        st.error(result_to_display)
    else:
        st.markdown(result_to_display)
        # Phần lưu vào bộ sưu tập chỉ hiện ra khi có kết quả thành công
        identified_name_from_result = result_to_display.split("**")[1] if "**" in result_to_display else None
        if is_logged_in and identified_name_from_result:
            st.markdown("---")
            st.subheader(f"Lưu '{identified_name_from_result}' vào bộ sưu tập")
            collections = st.session_state.get("collections", {})
            if not collections:
                st.info("Bạn chưa có bộ sưu tập nào. Hãy tạo ở thanh công cụ bên trái.")
            else:
                col1, col2 = st.columns([2,1])
                with col1:
                    selected_collection = st.selectbox("Chọn bộ sưu tập:", options=list(collections.keys()), key="collection_selector")
                with col2:
                    st.write("") 
                    st.write("")
                    if st.button("Thêm thuốc", use_container_width=True):
                        user_info = st.session_state.user_info
                        message = utils.add_drug_to_collection(firebase_db, user_info, selected_collection, identified_name_from_result)
                        st.toast(message)
                        _, collections_new = utils.load_user_data(firebase_db, user_info)
                        st.session_state.collections = collections_new
                        st.rerun()


# --- SIDEBAR ---
with st.sidebar:
    st.header("Lịch sử tra cứu")
    
    # Khung chứa có thanh cuộn áp dụng cho cả khách và người dùng
    history_container = st.container(height=300)
    with history_container:
        if not st.session_state.history:
            st.info("Chưa có thuốc nào được tra cứu.")
        else:
            for drug in st.session_state.history:
                col1, col2 = st.columns([0.8, 0.2])
                with col1:
                    if st.button(drug, key=f"history_{drug}", use_container_width=True):
                        run_lookup(drug)
                with col2:
                    if is_logged_in:
                        with st.popover("➕", use_container=True):
                            collections = st.session_state.get("collections", {})
                            if not collections:
                                st.write("Chưa có bộ sưu tập.")
                            else:
                                for coll_name in collections.keys():
                                    if st.button(f"Thêm vào '{coll_name}'", key=f"add_{drug}_to_{coll_name}"):
                                        user_info = st.session_state.user_info
                                        message = utils.add_drug_to_collection(firebase_db, user_info, coll_name, drug)
                                        st.toast(message)
                                        _, collections_new = utils.load_user_data(firebase_db, user_info)
                                        st.session_state.collections = collections_new
                                        st.rerun()

    st.markdown("---")

    # --- PHẦN BỘ SƯU TẬP ---
    if is_logged_in:
        st.header("Bộ sưu tập")
        def handle_create_collection():
            coll_name = st.session_state.new_collection_input
            success, message = utils.create_new_collection(firebase_db, st.session_state.user_info, coll_name)
            if success:
                st.success(message)
                st.session_state.collections[coll_name] = []
                st.session_state.new_collection_input = ""
            else:
                st.error(message)

        st.text_input("Tên bộ sưu tập mới:", key="new_collection_input")
        st.button("Tạo mới", on_click=handle_create_collection)

        collections = st.session_state.get("collections", {})
        for name, drugs in collections.items():
            with st.expander(f"{name} ({len(drugs)}/{utils.DRUGS_PER_COLLECTION_LIMIT} thuốc)"):
                if not drugs:
                    st.write("Bộ sưu tập này trống.")
                else:
                    for drug in drugs:
                        if st.button(drug, key=f"collection_{name}_{drug}"):
                            run_lookup(drug)
        st.markdown(f"Đã tạo {len(collections)}/{utils.COLLECTION_LIMIT} bộ sưu tập.")
        st.markdown("---")

    with st.container(border=True):
        st.write("**Bạn có ý tưởng để cải thiện ứng dụng?**")
        st.link_button("Gửi phản hồi ngay!", url="https://forms.gle/M44GDS4hJ7LpY7b98")

    if is_logged_in:
        st.header("Truy cập Pro")
        if st.session_state.get("pro_access"):
            st.success("Bạn đã có quyền truy cập Pro.")
        else:
            pro_code_input = st.text_input("Nhập mã truy cập Pro:", type="password")
            if st.button("Xác thực"):
                is_valid, message = utils.verify_code(pro_code_input)
                if is_valid:
                    st.success(message)
                    st.rerun()
                else:
                    st.error(message)
