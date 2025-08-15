import streamlit as st
import google.generativeai as genai

# Import các file chúng ta đã tách
import auth
import utils

# --- KHỞI TẠO CÁC DỊCH VỤ ---
# Khởi tạo Firebase App và các dịch vụ con
firebase_app = auth.initialize_firebase_app()
if not firebase_app:
    st.stop() # Dừng ứng dụng nếu không kết nối được Firebase

firebase_auth = firebase_app.auth()
firebase_db = firebase_app.database() # Lấy đối tượng database để tương tác

# Cấu hình Google AI
try:
    genai.configure(api_key=st.secrets["GOOGLE_API_KEY"])
except (FileNotFoundError, KeyError):
    st.error("LỖI: Vui lòng cấu hình GOOGLE_API_KEY trong secrets.toml.")
    st.stop()

# --- KIỂM TRA TRẠNG THÁI BẢO TRÌ ---
is_maintenance = st.secrets.get("maintenance_mode", False)
if is_maintenance:
    st.set_page_config(page_title="Bảo trì", page_icon="🛠️")
    st.title("🛠️ Dược Điển AI đang được bảo trì")
    message = st.secrets.get("maintenance_message", "Ứng dụng đang được cập nhật. Vui lòng quay lại sau.")
    st.info(message)
    st.stop()

# --- KHỞI TẠO TRẠNG THÁI PHIÊN ---
if 'history' not in st.session_state: st.session_state.history = []
if 'pro_access' not in st.session_state: st.session_state.pro_access = False
if 'user_data_loaded' not in st.session_state: st.session_state.user_data_loaded = False
if 'collections' not in st.session_state: st.session_state.collections = {}
if 'last_drug_searched' not in st.session_state: st.session_state.last_drug_searched = None

# --- HÀM LOGIC TRUNG TÂM ---
def run_lookup(drug_name):
    st.session_state.last_drug_searched = None # Reset thuốc vừa tra cứu
    try:
        is_pro = st.session_state.get("pro_access", False)
        final_result = utils.get_drug_info(drug_name, is_pro_user=is_pro)
        
        if not final_result.startswith("❌ Lỗi:"):
            st.markdown(final_result)
            user_info = st.session_state.get("user_info")
            hoat_chat_da_nhan_dien = final_result.split("**")[1]
            st.session_state.last_drug_searched = hoat_chat_da_nhan_dien # Lưu lại thuốc vừa tra cứu thành công
            
            if user_info: # Nếu người dùng đã đăng nhập
                utils.save_drug_to_history(firebase_db, user_info, hoat_chat_da_nhan_dien)
                st.session_state.history = utils.load_user_history(firebase_db, user_info)
            else: # Nếu là khách
                if drug_name not in st.session_state.history:
                    st.session_state.history.insert(0, drug_name)
                    if len(st.session_state.history) > 10:
                        st.session_state.history.pop()
        else:
            st.error(final_result)
    except Exception as e:
        st.error("💥 Lỗi không xác định.")
        st.exception(e)

# --- BẮT ĐẦU GIAO DIỆN ---
st.set_page_config(page_title="Dược Điển AI", page_icon="💊")

# Hiển thị form đăng nhập và lấy trạng thái
is_logged_in = auth.display_auth_forms(firebase_auth)

# --- TẢI DỮ LIỆU NGƯỜI DÙNG KHI ĐĂNG NHẬP ---
if is_logged_in and not st.session_state.user_data_loaded:
    user_info = st.session_state.user_info
    with st.spinner("Đang tải dữ liệu của bạn..."):
        st.session_state.history = utils.load_user_history(firebase_db, user_info)
        st.session_state.collections = utils.load_user_collections(firebase_db, user_info)
    st.session_state.user_data_loaded = True

# --- GIAO DIỆN CHÍNH ---
st.title("Dược Điển AI 💊")
st.caption("Dự án được phát triển bởi group CÂCK và AI")

# --- KHUNG NHẬP LIỆU CHÍNH ---
drug_name_input = st.text_input("Nhập tên thuốc (biệt dược hoặc hoạt chất):", key="main_input")
if st.button("Tra cứu"):
    if not drug_name_input:
        st.warning("Vui lòng nhập tên thuốc trước khi tra cứu.")
    else:
        run_lookup(drug_name_input)

# --- KHU VỰC LƯU VÀO BỘ SƯU TẬP (CHỈ HIỆN KHI CẦN) ---
if is_logged_in and st.session_state.last_drug_searched:
    st.markdown("---")
    st.subheader(f"Lưu '{st.session_state.last_drug_searched}' vào bộ sưu tập")
    
    collections = st.session_state.get("collections", {})
    if not collections:
        st.info("Bạn chưa có bộ sưu tập nào. Hãy tạo ở thanh công cụ bên trái.")
    else:
        # --- LOGIC CHO NÚT "THÊM THUỐC" DÙNG CALLBACK ---
        def handle_add_drug_to_collection():
            user_info = st.session_state.user_info
            drug_to_add = st.session_state.last_drug_searched
            # Lấy collection được chọn từ session_state
            selected_collection = st.session_state.get("collection_selector")

            if utils.add_drug_to_collection(firebase_db, user_info, selected_collection, drug_to_add):
                st.success(f"Đã thêm '{drug_to_add}' vào '{selected_collection}'.")
                st.session_state.collections = utils.load_user_collections(firebase_db, user_info)
            else:
                st.warning(f"'{drug_to_add}' đã có trong '{selected_collection}'.")

        col1, col2 = st.columns([2,1])
        with col1:
            st.selectbox("Chọn bộ sưu tập:", options=list(collections.keys()), key="collection_selector")
        with col2:
            st.write("") 
            st.write("")
            st.button("Thêm thuốc", on_click=handle_add_drug_to_collection, use_container_width=True)

# --- SIDEBAR ---
with st.sidebar:
    st.header("Lịch sử tra cứu")
    if not st.session_state.history:
        st.info("Chưa có thuốc nào được tra cứu.")
    else:
        for drug in st.session_state.history:
            if st.button(drug, key=f"history_{drug}", use_container_width=True):
                st.session_state.main_input = drug
                st.rerun()

    st.markdown("---")

    # --- PHẦN BỘ SƯU TẬP TRÊN SIDEBAR ---
    if is_logged_in:
        st.header("Bộ sưu tập")
        
        # --- LOGIC "TẠO MỚI" DÙNG CALLBACK ---
        def handle_create_collection():
            collection_name_to_create = st.session_state.new_collection_input
            if not collection_name_to_create or collection_name_to_create.isspace():
                st.error("Tên bộ sưu tập không được để trống.")
                return

            user_info = st.session_state.user_info
            success, message = utils.create_new_collection(firebase_db, user_info, collection_name_to_create)
            if success:
                st.success(message)
                st.session_state.collections = utils.load_user_collections(firebase_db, user_info)
                # Xóa chữ trong ô input sau khi tạo thành công - an toàn trong callback
                st.session_state.new_collection_input = ""
            else:
                st.error(message)

        st.text_input("Tên bộ sưu tập mới:", key="new_collection_input")
        st.button("Tạo mới", on_click=handle_create_collection)

        collections = st.session_state.get("collections", {})
        if not collections:
            st.info("Chưa có bộ sưu tập nào.")
        else:
            for name, drugs in collections.items():
                with st.expander(f"{name} ({len(drugs)} thuốc)"):
                    if not drugs:
                        st.write("Bộ sưu tập này trống.")
                    else:
                        for drug in drugs:
                            if st.button(drug, key=f"collection_{name}_{drug}", use_container_width=True):
                                st.session_state.main_input = drug
                                st.rerun()
        st.markdown("---")

    with st.container(border=True):
        st.write("**Bạn có ý tưởng để cải thiện ứng dụng?**")
        st.link_button("Gửi phản hồi ngay!", url="https://forms.gle/M44GDS4hJ7LpY7b98", help="Mở form góp ý trong một tab mới")

    if is_logged_in:
        st.header("Truy cập Pro")
        if st.session_state.get("pro_access"):
            st.success("Bạn đã có quyền truy cập Pro.")
        else:
            pro_code_input = st.text_input("Nhập mã truy cập Pro:", type="password", help="Nhập mã của bạn và bấm nút Xác thực.")
            if st.button("Xác thực"):
                is_valid, message = utils.verify_code(pro_code_input)
                if is_valid:
                    st.success(message)
                    st.rerun()
                else:
                    st.error(message)
