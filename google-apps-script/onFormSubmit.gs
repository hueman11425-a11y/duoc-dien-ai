function onFormSubmit(e) {
  const recipientEmail = 'abc@gmail.com'; // Thay gmail phù hợp
  let isTest = false;

  // Nếu chạy trực tiếp → dùng dữ liệu giả lập
  if (!e || !e.values) {
    Logger.log("⚠️ Không có e.values → đang chạy TEST.");
    isTest = true;
    e = {
      values: [
        "15/08/2025 18:42:27", // Cột A - Dấu thời gian
        "Báo lỗi ứng dụng",     // Cột B - Loại phản hồi
        "Ứng dụng bị treo khi lưu", // Cột C - Nội dung chi tiết
        "https://drive.google.com/file/d/FAKE_FILE_ID/view?usp=drive_link", // Cột D - Hình ảnh báo lỗi
        "abc@example.com"       // Cột E - Thông tin liên hệ
      ]
    };
  }

  const row = e.values;

  // Mapping cột (bắt đầu từ 0)
  const timestamp       = row[0] || "";
  const feedbackType    = row[1] || "Không có";
  const feedbackDetails = row[2] || "Không có";
  const fileLink        = row[3] || "";
  const contactInfo     = row[4] || "Ẩn danh";

  // Soạn tiêu đề email
  const emailSubject = `${isTest ? "[TEST] " : ""}🔔 Phản hồi mới cho Dược Điển AI: ${feedbackType}`;

  // Soạn nội dung email (HTML)
  let emailBody = `
    <div style="font-family:Segoe UI, sans-serif; font-size:14px; color:#222;">
      <h2 style="color:#0078D4;">🔔 Báo cáo phản hồi mới</h2>
      <table style="border-collapse:collapse;">
        <tr><td><b>Dấu thời gian:</b></td><td>${timestamp}</td></tr>
        <tr><td><b>Loại phản hồi:</b></td><td>${feedbackType}</td></tr>
        <tr><td valign="top"><b>Nội dung chi tiết:</b></td>
            <td>${feedbackDetails.replace(/\n/g, '<br>')}</td></tr>
        <tr><td><b>Liên hệ:</b></td><td>${contactInfo}</td></tr>
      </table>
  `;

  if (fileLink) {
    emailBody += `
      <p><b>📎 Tệp đính kèm:</b> <a href="${fileLink}" target="_blank">Xem tệp</a></p>
    `;
  }

  emailBody += `<p style="margin-top:20px;">Chúc bạn một ngày làm việc hiệu quả!<br>— Bot Báo cáo Dược Điển AI</p></div>`;

  // Gửi email
  MailApp.sendEmail({
    to: recipientEmail,
    subject: emailSubject,
    htmlBody: emailBody,
    name: 'Bot Báo cáo Dược Điển AI'
  });

  Logger.log("✅ Email đã được gửi.");
}
