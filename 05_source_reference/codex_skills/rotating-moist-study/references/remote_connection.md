# 远端SSH连接 / Remote SSH connection

## Verified active profile

Verified replacement platform on 2026-08-15:

```text
Host xh5
HostName xh5.hpccube.com
User jiasenzhang
Port 65061
```

The 14-case continuation bundle is located at:

```text
/work/home/jiasenzhang/rotating_moist_migration_bundle_20260815
```

This platform uses Slurm partition `xhacnormalc`. Prepare and validate jobs
with `04_platform_scripts/submit_all.sh`; do not add `--submit` without the
user's explicit approval.

Verified on 2026-08-05 from the current Windows workstation:

```text
Host c01n0011
HostName phssh.hnaicc.cn
User shu_zhangjs
Port 13470
ForwardAgent yes
IdentityFile ~/.ssh/74620990_rsa
```

Successful verification returned:

```text
hostname: c01n0011
user: shu_zhangjs
home: /share/org/SHUTUANL/shu_zhangjs
```

The private key is a credential. Never place it inside a skill, transfer ZIP,
Git repository, chat message, plot folder, or remote postprocessing archive.
Prefer issuing a separate key for each physical device through the cluster
portal. If the cluster explicitly permits reuse of the existing key, the user
must transfer it separately through a secure channel.

## Configure another Windows Codex device

1. Install or enable Windows OpenSSH Client.
2. Put the authorized private key at:

```text
C:\Users\<WindowsUser>\.ssh\74265428_rsa
```

3. Create or edit:

```text
C:\Users\<WindowsUser>\.ssh\config
```

4. Add the verified profile above. The portable `~/.ssh/...` form avoids
hard-coding the old computer's Windows username.
5. Restrict the private-key ACL if OpenSSH reports bad permissions:

```powershell
icacls "$HOME\.ssh\74265428_rsa" /inheritance:r
icacls "$HOME\.ssh\74265428_rsa" /grant:r "$env:USERNAME:(R)"
```

6. Validate configuration without connecting:

```powershell
ssh -G c01n0006 | Select-String "hostname|user|port|identityfile"
```

7. Test the real connection:

```powershell
ssh c01n0011 "hostname; whoami; pwd"
```

Expected output contains `c01n0011`, `shu_zhangjs`, and
`/share/org/SHUTUANL/shu_zhangjs`.

## Codex operating rules after connection

- Use `ssh c01n0011` or `scp c01n0011:<remote-path> <local-path>`; do not
replace the alias with an old node name or port from memory.
- Process large HDF5 movie fields remotely. Download only reduced CSV, NPY,
NPZ, JSON, PNG, PDF, XMF, or compact HDF5 products.
- Quote remote paths containing spaces, especially `rainy model`.
- Do not submit jobs unless the user explicitly approves the prepared cases.
- When reading running cases, record the scan time because movie and online
diagnostic files continue to advance.

Main study root:

```text
/share/org/SHUTUANL/shu_zhangjs/rainy model/ns/transition_study/beta1
```

AR10 study root:

```text
/share/org/SHUTUANL/shu_zhangjs/rainy model/ns/transition_study/beta1/lowrestest/aspect_ratio_study
```

## Troubleshooting

- `Permission denied (publickey)`: verify username, identity path, key ACL,
  key expiration, and whether the key is authorized for this device/account.
- Timeout or connection refused: verify `phssh.hnaicc.cn:13366`, network/VPN,
  and the cluster portal status.
- Host-key warning: verify the new fingerprint through the cluster portal or
  administrator. Do not blindly delete `known_hosts` entries.
- Codex can run `ssh` only if the Codex process has access to the current
  Windows user's `.ssh` directory and the workspace policy permits network
  commands.
