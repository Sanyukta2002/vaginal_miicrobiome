#!/usr/bin/env bash -C -e -u -o pipefail
printf "%s %s\n" 1/S2007001548_1.fastq.gz S2007001548_1.gz 2/S2007001548_2.fastq.gz S2007001548_2.gz | while read old_name new_name; do
    [ -f "${new_name}" ] || ln -s $old_name $new_name
done

fastqc \
    --quiet \
    --threads 2 \
    --memory 1536 \
    S2007001548_1.gz S2007001548_2.gz

# capture process environment
set +u
set +e
cd "$NXF_TASK_WORKDIR"

nxf_eval_cmd() {
    {
        IFS=$'\n' read -r -d '' "${1}";
        IFS=$'\n' read -r -d '' "${2}";
        (IFS=$'\n' read -r -d '' _ERRNO_; return ${_ERRNO_});
    } < <((printf '\0%s\0%d\0' "$(((({ shift 2; "${@}"; echo "${?}" 1>&3-; } | tr -d '\0' 1>&4-) 4>&2- 2>&1- | tr -d '\0' 1>&4-) 3>&1- | exit "$(cat)") 4>&1-)" "${?}" 1>&2) 2>&1)
}

echo '' > .command.env
#
nxf_eval_cmd STDOUT STDERR bash -c "fastqc --version | sed \"/FastQC v/!d; s/.*v//\""
status=$?
if [ $status -eq 0 ]; then
  echo nxf_out_eval_1="$STDOUT" >> .command.env
  echo /nxf_out_eval_1/=exit:0 >> .command.env
else
  echo nxf_out_eval_1="$STDERR" >> .command.env
  echo /nxf_out_eval_1/=exit:$status >> .command.env
fi
