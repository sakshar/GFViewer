(function () {
  "use strict";
  var form = document.getElementById("gfv-form");
  if (!form) return;
  var msg = document.getElementById("gfv-form-msg");
  var btn = form.querySelector("button[type=submit]");

  form.addEventListener("submit", function (ev) {
    ev.preventDefault();
    msg.textContent = "";
    msg.className = "msg";
    btn.disabled = true;
    btn.textContent = "Uploading…";

    fetch(form.action, { method: "POST", body: new FormData(form) })
      .then(function (r) {
        return r.json().then(function (body) { return { ok: r.ok, body: body }; });
      })
      .then(function (res) {
        if (!res.ok) {
          var text = res.body.error || "Submission failed.";
          if (res.body.hints && res.body.hints.length) {
            text += " " + res.body.hints.join(" ");
          }
          throw new Error(text);
        }
        window.location.href = "/jobs/" + res.body.job_id;
      })
      .catch(function (err) {
        msg.textContent = err.message;
        msg.className = "msg error";
        btn.disabled = false;
        btn.textContent = "Generate figure";
      });
  });
})();
