import streamlit as st
from utils import drugs, firebase, constants, prescription
import markdown2
import win32clipboard
import re

# ===========================
# HÀM HỖ TRỢ: Chuẩn hóa danh sách số và Markdown → HTML
# ===========================
def normalize_numbered_list(text: str) -> str:
    """
    Chuẩn hóa các dòng bắt đầu bằng số + '.'.
    Bỏ các dòng trống.
    """
    lines = text.splitlines()
    new_lines = []
    for line in lines:
        line_strip = line.strip()
        if not line_strip:
            continue
        match = re.match(r'^(\d+)\.\s*(\S.+)$', line_strip)
        if match:
            new_lines.append(f"{match.group(1)}. {match.group(2)}")
        else:
            new_lines.append(line_strip)
    return "\n".join(new_lines)

def generate_styled_html(markdown_text: str) -> str:
    """
    Convert Markdown → HTML đẹp, giữ danh sách số,
    loại bỏ dòng trống, áp dụng font Segoe UI inline.
    """
    lines = [l for l in markdown_text.splitlines() if l.strip()]
    markdown_text = "\n".join(lines)
    markdown_text = normalize_numbered_list(markdown_text)
    html_body = markdown2.markdown(
        markdown_text,
        extras=["fenced-code-blocks", "tables", "break-on-newline"]
    )
    html = f'<div style="font-family:\'Segoe UI\'; font-size:12pt; color:#222; line-height:1.5;">{html_body}</div>'
    return html

def copy_html_to_clipboard(html: str):
    """
    Copy HTML vào clipboard (Windows) chuẩn CF_HTML
    """
    html_bytes = html.encode("utf-8")
    win32clipboard.OpenClipboard()
    win32clipboard.EmptyClipboard()
    fmt = win32clipboard.RegisterClipboardFormat("HTML Format")
    win32clipboard.SetClipboardData(fmt, html_bytes)
    win32clipboard.CloseClipboard()

# ===========================
# TRANG TRA CỨU THUỐC
# ===========================
def render_lookup_page():
    st.header("Tra cứu Dược điển 💊")
    is_logged_in = st.session_state.user_info is not None
    firebase_db = st.session_state.firebase_db
    user_info = st.session_state.user_info

    def run_lookup(drug_name: str):
        st.session_state.query_result = None
        is_pro = st.session_state.get("pro_access", False)
        with st.spinner(f"Đang tra cứu '{drug_name}'..."):
            api_result, identified_name = drugs.get_drug_info_from_api(drug_name, is_pro)
            st.session_state.query_result = api_result
            st.session_state.identified_name = identified_name  # ← lưu vào session_state

            if identified_name and "không được nhận dạng" not in api_result:
                if is_logged_in:
                    updated_history = firebase.save_new_result(firebase_db, user_info, identified_name, api_result)
                    if updated_history is not None:
                        st.session_state.history = updated_history
                else:
                    st.session_state.guest_cache[identified_name] = api_result
                    if identified_name not in st.session_state.history:
                        st.session_state.history.insert(0, identified_name)
                        if len(st.session_state.history) > 10:
                            drug_to_remove = st.session_state.history.pop()
                            st.session_state.guest_cache.pop(drug_to_remove, None)

    search_query = st.text_input("Nhập tên thuốc:", key="lookup_input_field")

    if st.button("Tra cứu"):
        if search_query:
            run_lookup(search_query)
        else:
            st.warning("Vui lòng nhập tên thuốc trước khi tra cứu.")

    if st.session_state.query_result:
        result_text = st.session_state.query_result
        if result_text.startswith("❌ Lỗi:"):
            st.error(result_text)
        else:
            lines = result_text.split("\n", 1)
            first_line = lines[0] if lines else ""
            rest = lines[1] if len(lines) > 1 else ""
            st.markdown(
                f"<div style=\"font-family:'Segoe UI'; font-size:18pt;font-weight:bold;color:#111\">{first_line}</div>",
                unsafe_allow_html=True
            )
            st.markdown(rest.replace("\n", "  \n"))

            if st.button("📋 Sao chép kết quả", key="copy_button"):
                try:
                    html_content = generate_styled_html(result_text)
                    copy_html_to_clipboard(html_content)
                    st.success("✅ Đã sao chép vào clipboard")
                except Exception as e:
                    st.error(f"❌ Không thể sao chép: {e}")

            # ===== Thêm thuốc vào bộ sưu tập =====
            if is_logged_in and st.session_state.get("identified_name"):
                drug_name = st.session_state.identified_name
                st.markdown("---")
                st.subheader(f"Lưu '{drug_name}' vào bộ sưu tập")
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
                            try:
                                msg = firebase.add_drug_to_collection(firebase_db, user_info, sel_coll, drug_name)
                                if msg is None:
                                    msg = f"Đã thêm '{drug_name}' vào '{sel_coll}'"
                                st.toast(msg)
                                _, collections_new = firebase.load_user_data(firebase_db, user_info)
                                st.session_state.collections = collections_new
                                st.success(f"✅ '{drug_name}' đã được thêm vào '{sel_coll}'")
                                st.rerun()
                            except Exception as e:
                                st.error(f"❌ Không thể thêm thuốc: {e}")

# ===========================
# TRANG PHÂN TÍCH ĐƠN THUỐC
# ===========================
def render_prescription_analysis_page():
    st.header("Phân tích Đơn thuốc 🩺")
    is_logged_in = st.session_state.user_info is not None

    if not is_logged_in:
        st.warning("Vui lòng đăng nhập để sử dụng tính năng này.")
        return

    st.subheader("Phần 1: Bối cảnh Bệnh nhân (Không bắt buộc)")
    col1, col2 = st.columns(2)
    with col1:
        conditions = st.text_area("Tình trạng bệnh lý nền:", height=150)
    with col2:
        allergies = st.text_area("Dị ứng thuốc đã biết:", height=150)

    st.subheader("Phần 2: Thông tin Đơn thuốc (Bắt buộc)")
    prescription_text = st.text_area("Dán nội dung đơn thuốc:", height=250)

    if st.button("Phân tích Đơn thuốc", type="primary"):
        if not prescription_text.strip():
            st.error("Vui lòng nhập thông tin đơn thuốc.")
        else:
            with st.spinner("AI đang phân tích..."):
                patient_context = f"- Bệnh lý nền: {conditions or 'Không có'}\n- Dị ứng: {allergies or 'Không có'}"
                result = prescription.get_prescription_analysis(
                    st.session_state.firebase_db, st.session_state.user_info,
                    patient_context, prescription_text
                )
                st.session_state.analysis_result = result

    if st.session_state.analysis_result:
        st.markdown("---")
        st.subheader("Kết quả Phân tích")
        result_text = st.session_state.analysis_result
        if result_text.startswith("❌"):
            st.error(result_text)
        else:
            st.markdown(result_text.replace("\n", "  \n"))

            if st.button("📋 Sao chép kết quả phân tích", key="copy_analysis_button"):
                try:
                    html_content = generate_styled_html(result_text)
                    copy_html_to_clipboard(html_content)
                    st.success("✅ Đã sao chép vào clipboard")
                except Exception as e:
                    st.error(f"❌ Không thể sao chép: {e}")
