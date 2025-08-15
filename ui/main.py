import streamlit as st
import markdown2
import win32clipboard

def generate_styled_html(markdown_text: str):
    # Tách dòng đầu tiên (tiêu đề)
    lines = markdown_text.split("\n", 1)
    first_line = lines[0] if lines else ""
    rest = lines[1] if len(lines) > 1 else ""

    # Chuyển phần còn lại sang HTML
    rest_html = markdown2.markdown(rest)

    # CSS: tiêu đề giữ Segoe UI, phần còn lại font mặc định
    css = """
    <style>
        body { font-family: Arial, sans-serif; font-size: 14pt; color: #222; }
        .drug-title { font-family: 'Segoe UI', sans-serif; font-size: 18pt; font-weight:bold; color:#111; margin-bottom:8px; }
        h1, h2, h3 { color: #111; }
        strong { font-weight: bold; }
        ul, ol { margin-left: 20px; }
        li { margin-bottom: 4px; }
        hr { border: 1px solid #ccc; }
    </style>
    """

    html = f"""
    <html>
    <head>{css}</head>
    <body>
        <div class="drug-title">{first_line}</div>
        {rest_html}
    </body>
    </html>
    """
    return html

def copy_html_to_clipboard(html: str):
    html_bytes = html.encode("utf-8")
    start_html = 0
    end_html = len(html_bytes)
    header = f"Version:0.9\r\nStartHTML:{start_html:010}\r\nEndHTML:{end_html:010}\r\nStartFragment:{start_html:010}\r\nEndFragment:{end_html:010}\r\n"
    full_content = header.encode("utf-8") + html_bytes

    win32clipboard.OpenClipboard()
    win32clipboard.EmptyClipboard()
    fmt = win32clipboard.RegisterClipboardFormat("HTML Format")
    win32clipboard.SetClipboardData(fmt, full_content)
    win32clipboard.CloseClipboard()

def render_main_interface(get_drug_info):
    drug_name = st.text_input("Nhập tên thuốc:")

    if st.button("Tra cứu"):
        if not drug_name:
            st.warning("Vui lòng nhập tên thuốc.")
        else:
            is_pro = st.session_state.get("pro_access", False)
            result = get_drug_info(drug_name, is_pro_user=is_pro)
            st.session_state["result_text"] = result
            st.session_state["last_query"] = drug_name
            if not result.startswith("❌ Lỗi:"):
                _update_history(drug_name)

    result = st.session_state.get("result_text")
    if result:
        # Tách dòng đầu tiên (tiêu đề thuốc) và phần còn lại
        lines = result.split("\n", 1)
        first_line = lines[0] if lines else ""
        rest = lines[1] if len(lines) > 1 else ""

        # Hiển thị tiêu đề với font Segoe UI
        st.markdown(f"""
        <div style="font-family: 'Segoe UI', sans-serif; font-size: 18pt; font-weight:bold; color:#111;">
            {first_line}
        </div>
        """, unsafe_allow_html=True)

        # Hiển thị phần còn lại bình thường
        st.markdown(rest)

        if st.button("📋 Sao chép kết quả", key="copy_button"):
            try:
                html_content = generate_styled_html(result)
                copy_html_to_clipboard(html_content)
                st.success("✅ Đã sao chép vào clipboard")
            except Exception as e:
                st.error(f"❌ Không thể sao chép: {e}")

def _update_history(drug_name: str):
    if "history" not in st.session_state:
        st.session_state["history"] = []
    if drug_name not in st.session_state.history:
        st.session_state.history.insert(0, drug_name)
        if len(st.session_state.history) > 10:
            st.session_state.history.pop()
