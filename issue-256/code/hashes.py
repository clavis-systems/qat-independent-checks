"""Write SHA256SUMS.txt for code/ and logs/ (and the circuit) so a write-up can cite exact files."""
import hashlib, os

root = "/w"
lines = []
for sub in ("code", "logs"):
    for name in sorted(os.listdir(os.path.join(root, sub))):
        p = os.path.join(root, sub, name)
        if os.path.isfile(p) and not name.endswith(".pyc"):
            lines.append(f"{hashlib.sha256(open(p, 'rb').read()).hexdigest()}  {sub}/{name}")
c = "/circ/49Q_OLE_circuit_L_3_b_0.25_delta0.15.qasm"
lines.append(f"{hashlib.sha256(open(c, 'rb').read()).hexdigest()}  (tracker) {os.path.basename(c)}")
open(os.path.join(root, "SHA256SUMS.txt"), "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
