template<class T>
static __device__ __host__ void PE_area2_grad(const T in1[], const T in2[],
                            const T in3[], T grad[])
{
  T t10;
  T t11;
  T t12;
  T t13;
  T t14;
  T t15;
  T t16;
  T t26;
  T t27;
  T t28;
  T t8;
  T t9;
  // PE_area_grad_only_func
  //     GRAD = PE_area_grad_only_func(IN1,IN2,IN3)
  //     This function was generated by the Symbolic Math Toolbox version 9.1.
  //     09-Mar-2023 18:46:33
  t8 = in2[0] + -in3[0];
  t9 = in2[1] + -in3[1];
  t10 = in2[2] + -in3[2];
  t11 = in2[0] + -in1[0];
  t12 = in2[1] + -in1[1];
  t13 = in2[2] + -in1[2];
  t14 = in3[0] + -in1[0];
  t15 = in3[1] + -in1[1];
  t16 = in3[2] + -in1[2];
  t26 = t11 * t15 + -(t12 * t14);
  t27 = t11 * t16 + -(t13 * t14);
  t28 = t12 * t16 + -(t13 * t15);
  grad[0] = t9 * t26 * 2.0F + t10 * t27 * 2.0F;
  grad[1] = t8 * t26 * -2.0F + t10 * t28 * 2.0F;
  grad[2] = t8 * t27 * -2.0F - t9 * t28 * 2.0F;
  grad[3] = t15 * t26 * 2.0F + t16 * t27 * 2.0F;
  grad[4] = t14 * t26 * -2.0F + t16 * t28 * 2.0F;
  grad[5] = t14 * t27 * -2.0F - t15 * t28 * 2.0F;
  grad[6] = t12 * t26 * -2.0F - t13 * t27 * 2.0F;
  grad[7] = t11 * t26 * 2.0F - t13 * t28 * 2.0F;
  grad[8] = t11 * t27 * 2.0F + t12 * t28 * 2.0F;
}