# 💊 Dược Điển AI

**Dược Điển AI** là một ứng dụng tra cứu thông tin dược lý thông minh, tích hợp AI (Google Generative AI) và dữ liệu PubMed, giúp người dùng nhanh chóng tra cứu:
- Thông tin hoạt chất và biệt dược
- Phân tích chuyên sâu cho người dùng Pro
- Tóm tắt các nghiên cứu lâm sàng mới nhất (trong 2 năm gần đây)

Ứng dụng được phát triển với [Streamlit](https://streamlit.io) và hỗ trợ kết nối Google Sheets để quản lý mã truy cập Pro.

---

## 🚀 Tính năng nổi bật
- 🔍 **Tra cứu thuốc**: Hỗ trợ nhập tên biệt dược hoặc hoạt chất.
- 🧠 **Nhận diện hoạt chất** bằng AI và phân tích thông tin dược lý.
- 📚 **Tích hợp PubMed**: Tự động lấy và tóm tắt các nghiên cứu mới nhất.
- 👑 **Chế độ Pro**: Truy cập phân tích chuyên sâu và dữ liệu nâng cao.
- 📝 **Lịch sử tra cứu**: Lưu lại 10 kết quả gần nhất.
- 📢 **Gửi phản hồi nhanh** ngay trên ứng dụng.

---

## 📂 Cấu trúc thư mục

```
duocdien_ai/
│
├── app.py                  # Điểm vào chính của Streamlit
├── core/                    # Logic xử lý và cấu hình
│   ├── config.py
│   ├── access_control.py
│   └── drug_lookup.py
├── ui/                      # Giao diện người dùng
│   ├── sidebar.py
│   └── main.py
├── prompts/                 # Các file prompt cho AI
│   ├── prompt_nhandien.txt
│   ├── prompt_regular.txt
│   ├── prompt_pro.txt
│   └── prompt_summary.txt
├── requirements.txt         # Danh sách thư viện cần thiết
└── README.md                # Tài liệu dự án
```

---

## 🔧 Cài đặt

1. **Clone dự án**
```bash
git clone https://github.com/<username>/duoc-dien-ai.git
cd duoc-dien-ai
```

2. **Tạo môi trường ảo và cài thư viện**
```bash
python -m venv venv
source venv/bin/activate   # macOS / Linux
venv\Scripts\activate      # Windows

pip install -r requirements.txt
```

3. **Cấu hình API key và thông tin bảo mật**
   - Tạo file `secrets.toml` trong thư mục `.streamlit/`:
```toml
# .streamlit/secrets.toml

GOOGLE_API_KEY = "Your_GG_API_KEY"
# --- CHẾ ĐỘ BẢO TRÌ ---
# Đặt là true để BẬT chế độ bảo trì, false để TẮT.
maintenance_mode = false

# Nội dung thông báo khi bảo trì
maintenance_message = "Dược Điển AI đang được nâng cấp lên phiên bản mới với nhiều tính năng hấp dẫn. Chúng tôi sẽ sớm quay trở lại. Xin cảm ơn!"


# --- CẤU HÌNH KẾT NỐI GOOGLE SHEETS ---
[connections.gsheets]
spreadsheet = "Your_Project_GGSheet_File_Name" # <--- ĐẢM BẢO ĐÂY LÀ TÊN CHÍNH XÁC CỦA FILE GOOGLE SHEETS BẠN ĐÃ TẠO

# Dán các giá trị từ file .json của bạn vào đây
[connections.gsheets.credentials]
type = "service_account" # <-- Chép từ file .json
project_id = "duoc-dien-ai" # <-- Chép từ file .json
private_key_id = "your_private_key_id" # <-- Chép từ file .json
private_key = "your_private_key" # <-- Chép toàn bộ chuỗi private key, bao gồm cả -----BEGIN PRIVATE KEY----- và -----END PRIVATE KEY-----\n
client_email = "your_client_email" # <-- Chép từ file .json
client_id = "your_client_id" # <-- Chép từ file .json
auth_uri = "https://accounts.google.com/o/oauth2/auth"
token_uri = "https://oauth2.googleapis.com/token"
auth_provider_x509_cert_url = "your_auth_provider_x509_cert_url"
client_x509_cert_url = "your_client_x509_cert_url" # <-- Chép từ file .json

# --- CẤU HÌNH MODEL AI ---
[models]
regular = "gemini-2.5-flash-lite"
pro = "gemini-2.5-flash-lite"

# --- API KEYS ---
[api_keys]
pubmed = "yours_api_keys"

# --- FIREBASE CONFIG ---
[firebase]
apiKey = "your_firebase_apiKey"
authDomain = "your_firebase_authDomain"
projectId = "your_firebase_projectId"
storageBucket = "your_firebase_storageBucket"
messagingSenderId = "your_firebase_messagingSenderId"
appId = "your_firebase_appId"
databaseURL = "your_firebase_databaseURL"
```

---

## ▶️ Chạy ứng dụng

```bash
streamlit run app.py
```

Truy cập **http://localhost:8501** để sử dụng ứng dụng.

---

## 📜 Giấy phép
Dự án được phát triển cho mục đích học tập và nghiên cứu. Không sử dụng thay thế cho tư vấn y khoa chuyên nghiệp.

---

## 🤝 Đóng góp
Nếu bạn muốn đóng góp hoặc báo lỗi:
- Tạo **issue** trên GitHub
- Gửi phản hồi trực tiếp từ ứng dụng
- Liên hệ nhóm phát triển

---
