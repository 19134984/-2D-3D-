# 解析 list-directed Fortran 输出；每条奇矩记录固定包含 13 个数值，
# 即使编译器把一条记录折成多行，也按数值字段重新拼接。
BEGIN {
    stride = 13
    token_count = 0
    record_count = 0
}

/^#/ { next }

{
    for (field = 1; field <= NF; field++) {
        value[token_count % stride] = $field + 0.0
        token_count++
        if (token_count % stride == 0) {
            time_tf = value[0]
            qx_rms = value[3]
            qy_rms = value[4]
            q_rms = value[5]
            q_max = value[6]

            if (record_count == 0) first_time = time_tf
            record_count++
            last_time = time_tf
            last_qx = qx_rms
            last_qy = qy_rms
            last_qrms = q_rms
            last_qmax = q_max

            if (q_max > max_qmax) {
                max_qmax = q_max
                max_qmax_time = time_tf
            }
            if (q_rms > max_qrms) {
                max_qrms = q_rms
                max_qrms_time = time_tf
            }

            if (time_tf >= 40.0) {
                tail40_count++
                tail40_qx_sum += qx_rms
                tail40_qy_sum += qy_rms
                tail40_qrms_sum += q_rms
                if (q_max > tail40_max_qmax) {
                    tail40_max_qmax = q_max
                    tail40_max_qmax_time = time_tf
                }
            }
            if (time_tf >= 250.0) {
                tail250_count++
                tail250_qx_sum += qx_rms
                tail250_qy_sum += qy_rms
                tail250_qrms_sum += q_rms
                if (q_max > tail250_max_qmax) {
                    tail250_max_qmax = q_max
                    tail250_max_qmax_time = time_tf
                }
            }
        }
    }
}

END {
    printf("records=%d first_tf=%.9g last_tf=%.9g\n", record_count, first_time, last_time)
    printf("global_max_qrms=%.16e at_tf=%.9g\n", max_qrms, max_qrms_time)
    printf("global_max_qmax=%.16e at_tf=%.9g\n", max_qmax, max_qmax_time)
    printf("last_qx=%.16e last_qy=%.16e last_qrms=%.16e last_qmax=%.16e\n", \
        last_qx, last_qy, last_qrms, last_qmax)
    if (tail40_count > 0) {
        printf("tf_ge_40_count=%d mean_qx=%.16e mean_qy=%.16e mean_qrms=%.16e max_qmax=%.16e at_tf=%.9g\n", \
            tail40_count, tail40_qx_sum/tail40_count, tail40_qy_sum/tail40_count, \
            tail40_qrms_sum/tail40_count, tail40_max_qmax, tail40_max_qmax_time)
    }
    if (tail250_count > 0) {
        printf("tf_ge_250_count=%d mean_qx=%.16e mean_qy=%.16e mean_qrms=%.16e max_qmax=%.16e at_tf=%.9g\n", \
            tail250_count, tail250_qx_sum/tail250_count, tail250_qy_sum/tail250_count, \
            tail250_qrms_sum/tail250_count, tail250_max_qmax, tail250_max_qmax_time)
    }
}
