"�k
VHostIDLE"IDLE(1�����B�@9�y���@A�����B�@I�y���@a��z!���?i��z!���?�Unknown
lHostConv2D"sequential/conv1d/conv1d(1B`��"�~@9B`��"�~@AB`��"�~@IB`��"�~@a� �,ǲ?i��>��I�?�Unknown
~HostReluGrad"(gradient_tape/sequential/conv1d/ReluGrad(1�O��nVm@9�O��nVm@A�O��nVm@I�O��nVm@a�+Sg�ˡ?ii��Mf�?�Unknown
�HostMaxPoolGrad":gradient_tape/sequential/max_pooling1d/MaxPool/MaxPoolGrad(1����xMm@9����xMm@A����xMm@I����xMm@a�1ơ?i�-�Ѱ��?�Unknown
yHostMatMul"%gradient_tape/sequential/dense/MatMul(1�t��k@9�t��k@A�t��k@I�t��k@a@�<s�?ig��ؐ�?�Unknown
�HostConv2DBackpropFilter";gradient_tape/sequential/conv1d/conv1d/Conv2DBackpropFilter(1�z�G=f@9�z�G=f@A�z�G=f@I�z�G=f@a�&���?iL�F�h�?�Unknown
uHostMaxPool" sequential/max_pooling1d/MaxPool(1���Q�d@9���Q�d@A���Q�d@I���Q�d@a�S���?i���
�0�?�Unknown
nHostBiasAdd"sequential/conv1d/BiasAdd(1y�&1�]@9y�&1�]@Ay�&1�]@Iy�&1�]@ah�W���?i���d��?�Unknown
{	HostMatMul"'gradient_tape/sequential/dense/MatMul_1(15^�I�T@95^�I�T@A5^�I�T@I5^�I�T@a�Fh>_�?iT�N��%�?�Unknown
�
HostResourceApplyAdam"$Adam/Adam/update_2/ResourceApplyAdam(1��|?5�P@9��|?5�P@A��|?5�P@I��|?5�P@a�ж�R��?i�����w�?�Unknown
oHost_FusedMatMul"sequential/dense/Relu(1�l���N@9�l���N@A�l���N@I�l���N@a����}3�?ikĜٰ��?�Unknown
�HostDataset"3Iterator::Model::ParallelMap::Zip[1]::ForeverRepeat(1)\����?@9)\����?@A%��C�<@I%��C�<@a��+6wq?ip�E���?�Unknown
�HostBiasAddGrad"3gradient_tape/sequential/conv1d/BiasAdd/BiasAddGrad(11�Z�;@91�Z�;@A1�Z�;@I1�Z�;@a[�^%�p?i{�b�u�?�Unknown
hHostRelu"sequential/conv1d/Relu(1     �7@9     �7@A     �7@I     �7@a@k��l?i���Pl"�?�Unknown
fHostGreaterEqual"GreaterEqual(1}?5^��1@9}?5^��1@A}?5^��1@I}?5^��1@a��4��me?iz%��7�?�Unknown
�HostResourceApplyAdam"$Adam/Adam/update_4/ResourceApplyAdam(1�&1��0@9�&1��0@A�&1��0@I�&1��0@a��T?�td?ibzP�NL�?�Unknown
�HostResourceApplyAdam"$Adam/Adam/update_5/ResourceApplyAdam(1���S�e/@9���S�e/@A���S�e/@I���S�e/@a�l6�c?i��'Z_�?�Unknown
dHostDataset"Iterator::Model(1�z�G�:@9�z�G�:@A��S�[/@I��S�[/@a�[c?i����_r�?�Unknown
lHostIteratorGetNext"IteratorGetNext(1=
ףp}*@9=
ףp}*@A=
ףp}*@I=
ףp}*@a����p`?i��-�p��?�Unknown
�HostResourceApplyAdam""Adam/Adam/update/ResourceApplyAdam(1V-���(@9V-���(@AV-���(@IV-���(@a[�y{��]?i�\kBj��?�Unknown
[HostAddV2"Adam/add(1��~j�4(@9��~j�4(@A��~j�4(@I��~j�4(@a��2�]]?i����?�Unknown
gHostStridedSlice"strided_slice(1P��n�&@9P��n�&@AP��n�&@IP��n�&@ab5mS�[?i�j.���?�Unknown
qHostDataset"Iterator::Model::ParallelMap(1���x�f&@9���x�f&@A���x�f&@I���x�f&@a��ԜI-[?i��8~��?�Unknown
^HostGatherV2"GatherV2(1������%@9������%@A������%@I������%@a��f�SrZ?im��b���?�Unknown
�HostDataset"=Iterator::Model::ParallelMap::Zip[0]::FlatMap[0]::Concatenate(1�O��n,@9�O��n,@A�V�#@I�V�#@a������W?ik`�����?�Unknown
tHost_FusedMatMul"sequential/dense_1/BiasAdd(1�"��~j"@9�"��~j"@A�"��~j"@I�"��~j"@a�N��_WV?i���j���?�Unknown
tHostAssignAddVariableOp"AssignAddVariableOp(1J+� @9J+�@AJ+� @IJ+�@a��Y�`�S?i��Ț���?�Unknown
�HostResourceApplyAdam"$Adam/Adam/update_3/ResourceApplyAdam(1H�z��@9H�z��@AH�z��@IH�z��@aB���FS?iPj8��?�Unknown
oHostReadVariableOp"Adam/ReadVariableOp(1�(\��u@9�(\��u@A�(\��u@I�(\��u@a�,1�1S?i�����?�Unknown
{HostMatMul"'gradient_tape/sequential/dense_1/MatMul(1V-�@9V-�@AV-�@IV-�@a ��FR?i�S[H��?�Unknown
`HostGatherV2"
GatherV2_1(1;�O���@9;�O���@A;�O���@I;�O���@a��;أ�Q?i�qG���?�Unknown
| HostSelect"(binary_crossentropy/logistic_loss/Select(1#��~j�@9#��~j�@A#��~j�@I#��~j�@a���`4nQ?i7�w�A�?�Unknown
v!HostDataset"!Iterator::Model::ParallelMap::Zip(1�Zd;�K@9�Zd;�K@AZd;�O@IZd;�O@a� _�hP?i�^'v�?�Unknown
�"HostReadVariableOp"'sequential/dense/BiasAdd/ReadVariableOp(1���x�@9���x�@A���x�@I���x�@a~Yk�*P?iR][�'�?�Unknown
�#HostBiasAddGrad"4gradient_tape/sequential/dense_1/BiasAdd/BiasAddGrad(1�I+�@9�I+�@A�I+�@I�I+�@aA�M}�N?i�^�bI/�?�Unknown
v$HostMul"%binary_crossentropy/logistic_loss/mul(1����M�@9����M�@A����M�@I����M�@a�A�l0N?iƋ�i�6�?�Unknown
�%HostReadVariableOp"4sequential/conv1d/conv1d/ExpandDims_1/ReadVariableOp(1�(\���@9�(\���@A�(\���@I�(\���@a��smI?i����0=�?�Unknown
Y&HostPow"Adam/Pow(1-���F@9-���F@A-���F@I-���F@a��~�H?i8gWC�?�Unknown
~'HostSelect"*binary_crossentropy/logistic_loss/Select_1(1-���F@9-���F@A-���F@I-���F@a��[�ubG?i	O֥/I�?�Unknown
�(HostResourceApplyAdam"$Adam/Adam/update_1/ResourceApplyAdam(1V-2@9V-2@AV-2@IV-2@aE��Z�IG?i5-O�?�Unknown
�)HostDynamicStitch"/gradient_tape/binary_crossentropy/DynamicStitch(1��K7��@9��K7��@A��K7��@I��K7��@a���B��F?it��J�T�?�Unknown
V*HostMean"Mean(1���K�@9���K�@A���K�@I���K�@a^'�O�}E?i>���Z�?�Unknown
�+HostTile"6gradient_tape/binary_crossentropy/weighted_loss/Tile_1(1j�t��@9j�t��@Aj�t��@Ij�t��@a�S���XE?iS�g_�?�Unknown
\,HostGreater"Greater(1�Zd;�@9�Zd;�@A�Zd;�@I�Zd;�@a�΃K�wD?iG��d�?�Unknown
�-HostDataset"MIterator::Model::ParallelMap::Zip[0]::FlatMap[0]::Concatenate[0]::TensorSlice(1�n���@9�n���@A�n���@I�n���@a���_D?imx��i�?�Unknown
V.HostSum"Sum_2(1'1�Z@9'1�Z@A'1�Z@I'1�Z@aw���]�C?i���F�n�?�Unknown
~/HostMaximum")gradient_tape/binary_crossentropy/Maximum(1V-�@9V-�@AV-�@IV-�@aS�
��B?im�E�:s�?�Unknown
�0HostSelect"8gradient_tape/binary_crossentropy/logistic_loss/Select_3(1����K@9����K@A����K@I����K@a!r`B?i�r��w�?�Unknown
�1Host	ZerosLike":gradient_tape/binary_crossentropy/logistic_loss/zeros_like(1�A`��"@9�A`��"@A�A`��"@I�A`��"@aJ٩��GB?i'�8�d|�?�Unknown
z2HostLog1p"'binary_crossentropy/logistic_loss/Log1p(1`��"��@9`��"��@A`��"��@I`��"��@a��2T�.B?i�鍝���?�Unknown
}3HostMatMul")gradient_tape/sequential/dense_1/MatMul_1(1�V-@9�V-@A�V-@I�V-@a��l��A?ir!i?]��?�Unknown
]4HostCast"Adam/Cast_1(1���K7@9���K7@A���K7@I���K7@a�<$uA?iA&�����?�Unknown
j5HostMean"binary_crossentropy/Mean(1�rh��|	@9�rh��|	@A�rh��|	@I�rh��|	@ad�L��>?i�����?�Unknown
V6HostAddN"AddN(1�� �rh	@9�� �rh	@A�� �rh	@I�� �rh	@a�����>?i��m\��?�Unknown
�7HostDataset"?Iterator::Model::ParallelMap::Zip[1]::ForeverRepeat::FromTensor(1!�rh��@9!�rh���?A!�rh��@I!�rh���?aSmO��=>?i�e�&$��?�Unknown
v8HostExp"%binary_crossentropy/logistic_loss/Exp(1H�z�G@9H�z�G@AH�z�G@IH�z�G@a[�֌�t=?i� ��Ҙ�?�Unknown
j9HostReadVariableOp"ReadVariableOp(1�Q���@9�Q���@A�Q���@I�Q���@a��)�;?i~�|L��?�Unknown
�:HostReadVariableOp")sequential/dense_1/BiasAdd/ReadVariableOp(1�G�z�@9�G�z�@A�G�z�@I�G�z�@a�5�;?i�n���?�Unknown
V;HostCast"Cast(1��Q��@9��Q��@A��Q��@I��Q��@akƗ,�Q;?i��.'��?�Unknown
X<HostEqual"Equal(1��Q��@9��Q��@A��Q��@I��Q��@akƗ,�Q;?i˔�l���?�Unknown
�=HostDataset"-Iterator::Model::ParallelMap::Zip[0]::FlatMap(1V-�0@9V-�0@AH�z�G@IH�z�G@a?2f;?i�kY��?�Unknown
�>HostBiasAddGrad"2gradient_tape/sequential/dense/BiasAdd/BiasAddGrad(1�z�G�@9�z�G�@A�z�G�@I�z�G�@a�� ,�:?iж�C��?�Unknown
X?HostCast"Cast_2(1��(\��@9��(\��@A��(\��@I��(\��@a{�X�V/9?i�Aɩi��?�Unknown
�@HostAddV2"3gradient_tape/binary_crossentropy/logistic_loss/add(1q=
ףp@9q=
ףp@Aq=
ףp@Iq=
ףp@a�|F��8?itr(���?�Unknown
�AHost	ZerosLike"<gradient_tape/binary_crossentropy/logistic_loss/zeros_like_1(1q=
ףp@9q=
ףp@Aq=
ףp@Iq=
ףp@a�|F��8?i������?�Unknown
�BHostGreaterEqual".binary_crossentropy/logistic_loss/GreaterEqual(1��Q�@9��Q�@A��Q�@I��Q�@a�>���h8?i�������?�Unknown
oCHostSigmoid"sequential/dense_1/Sigmoid(1�n���@9�n���@A�n���@I�n���@a�@1�8?iHY���?�Unknown
vDHostNeg"%binary_crossentropy/logistic_loss/Neg(1��ʡE�@9��ʡE�@A��ʡE�@I��ʡE�@a�����7?i�L���?�Unknown
vEHostAssignAddVariableOp"AssignAddVariableOp_2(1��n��@9��n��@A��n��@I��n��@a2l�=��6?iw�:C��?�Unknown
�FHostBroadcastTo"-gradient_tape/binary_crossentropy/BroadcastTo(1�(\���@9�(\���@A�(\���@I�(\���@a!^y��5?i�i{8��?�Unknown
~GHostAssignAddVariableOp"Adam/Adam/AssignAddVariableOp(1o��ʡ@9o��ʡ@Ao��ʡ@Io��ʡ@a`p��c5?i�����?�Unknown
�HHostFloorDiv"*gradient_tape/binary_crossentropy/floordiv(1F����x@9F����x@AF����x@IF����x@a���&325?iM�K>���?�Unknown
}IHostReluGrad"'gradient_tape/sequential/dense/ReluGrad(1�O��n@9�O��n@A�O��n@I�O��n@a��-?��4?i�s�!��?�Unknown
�JHostSelect"8gradient_tape/binary_crossentropy/logistic_loss/Select_2(1P��n� @9P��n� @AP��n� @IP��n� @a�S�4?irhF����?�Unknown
tKHostReadVariableOp"Adam/Cast/ReadVariableOp(1��C�l��?9��C�l��?A��C�l��?I��C�l��?a�ͨ�#Z3?i���C��?�Unknown
xLHostCast"&gradient_tape/binary_crossentropy/Cast(1�$��C�?9�$��C�?A�$��C�?I�$��C�?a���d��2?i
Wm��?�Unknown
�MHost
Reciprocal":gradient_tape/binary_crossentropy/logistic_loss/Reciprocal(1��/�$�?9��/�$�?A��/�$�?I��/�$�?a�h���H2?i7�7���?�Unknown
�NHostSelect"6gradient_tape/binary_crossentropy/logistic_loss/Select(1m������?9m������?Am������?Im������?an���&2?i������?�Unknown
bOHostDivNoNan"div_no_nan_1(1+���?9+���?A+���?I+���?a���1?i�=�)��?�Unknown
wPHostReadVariableOp"div_no_nan_1/ReadVariableOp(1J+��?9J+��?AJ+��?IJ+��?a�&�tLn0?i�ۉ�6��?�Unknown
vQHostSub"%binary_crossentropy/logistic_loss/sub(1���S��?9���S��?A���S��?I���S��?a�b�b/?i�Ӫ�(��?�Unknown
�RHostSum"5gradient_tape/binary_crossentropy/logistic_loss/Sum_1(1����S�?9����S�?A����S�?I����S�?a񃣉 �.?i
n�N��?�Unknown
�SHostMul"3gradient_tape/binary_crossentropy/logistic_loss/mul(1㥛� ��?9㥛� ��?A㥛� ��?I㥛� ��?a%��}=�-?i�L�����?�Unknown
XTHostCast"Cast_3(1��ʡE�?9��ʡE�?A��ʡE�?I��ʡE�?a s��+?i�3+ѣ��?�Unknown
rUHostAdd"!binary_crossentropy/logistic_loss(1��ʡE�?9��ʡE�?A��ʡE�?I��ʡE�?a s��+?i�T��?�Unknown
[VHostPow"
Adam/Pow_1(1+�����?9+�����?A+�����?I+�����?a�_���*?iL�Ӈ���?�Unknown
XWHostCast"Cast_4(1�ʡE���?9�ʡE���?A�ʡE���?I�ʡE���?a-!�x��*?i�e[����?�Unknown
vXHostSum"%binary_crossentropy/weighted_loss/Sum(1� �rh��?9� �rh��?A� �rh��?I� �rh��?a�� %�'?ikg�r$��?�Unknown
�YHostReadVariableOp"(sequential/dense_1/MatMul/ReadVariableOp(1Zd;�O��?9Zd;�O��?AZd;�O��?IZd;�O��?aq��-�'?iE1~����?�Unknown
�ZHostReadVariableOp"(sequential/conv1d/BiasAdd/ReadVariableOp(1sh��|?�?9sh��|?�?Ash��|?�?Ish��|?�?a��>��Y'?i5շ���?�Unknown
�[HostReadVariableOp"&sequential/dense/MatMul/ReadVariableOp(1��� �r�?9��� �r�?A��� �r�?I��� �r�?a�˖�Oa&?i����{��?�Unknown
v\HostReadVariableOp"Adam/Cast_2/ReadVariableOp(1;�O��n�?9;�O��n�?A;�O��n�?I;�O��n�?aL��W\&?i+0l���?�Unknown
~]HostRealDiv")gradient_tape/binary_crossentropy/truediv(1Zd;�O�?9Zd;�O�?AZd;�O�?IZd;�O�?a)��c� %?icIVt1��?�Unknown
v^HostReadVariableOp"Adam/Cast_3/ReadVariableOp(1�x�&1�?9�x�&1�?A�x�&1�?I�x�&1�?a� L�r#?ic�h��?�Unknown
�_HostCast"3binary_crossentropy/weighted_loss/num_elements/Cast(1m������?9m������?Am������?Im������?an���&"?i�����?�Unknown
�`HostMul"5gradient_tape/binary_crossentropy/logistic_loss/mul_1(1w��/��?9w��/��?Aw��/��?Iw��/��?a�Ϯ�!?ir�7���?�Unknown
vaHostAssignAddVariableOp"AssignAddVariableOp_3(1�$��C�?9�$��C�?A�$��C�?I�$��C�?a�}'ࠉ ?i�t�Ѫ��?�Unknown
�bHostNeg"7gradient_tape/binary_crossentropy/logistic_loss/sub/Neg(1�x�&1�?9�x�&1�?A�x�&1�?I�x�&1�?aʱ^�'?ij����?�Unknown
TcHostMul"Mul(1bX9���?9bX9���?AbX9���?IbX9���?a�:n�?i��Y-q��?�Unknown
�dHostMul"7gradient_tape/binary_crossentropy/logistic_loss/mul/Mul(1�n����?9�n����?A�n����?I�n����?a�����o?i(��D��?�Unknown
veHostAssignAddVariableOp"AssignAddVariableOp_1(1�E�����?9�E�����?A�E�����?I�E�����?a	$�E?i������?�Unknown
�fHostSum"3gradient_tape/binary_crossentropy/logistic_loss/Sum(1P��n��?9P��n��?AP��n��?IP��n��?a#64�Q�?iS�=����?�Unknown
`gHostDivNoNan"
div_no_nan(1�Zd;��?9�Zd;��?A�Zd;��?I�Zd;��?aWxz��?i'/�Ȗ��?�Unknown
�hHostDivNoNan"@gradient_tape/binary_crossentropy/weighted_loss/value/div_no_nan(1bX9���?9bX9���?AbX9���?IbX9���?a�Gf@9�?iY2|*,��?�Unknown
�iHostSum"7gradient_tape/binary_crossentropy/logistic_loss/sub/Sum(1w��/��?9w��/��?Aw��/��?Iw��/��?a�Ϯ�?i֨!;���?�Unknown
�jHostNeg"3gradient_tape/binary_crossentropy/logistic_loss/Neg(1��Q���?9��Q���?A��Q���?I��Q���?a�p|a�?i��-�=��?�Unknown
}kHostDivNoNan"'binary_crossentropy/weighted_loss/value(1�z�G��?9�z�G��?A�z�G��?I�z�G��?a��<��M?i�;-6���?�Unknown
ylHostReadVariableOp"div_no_nan_1/ReadVariableOp_1(1�/�$�?9�/�$�?A�/�$�?I�/�$�?ad�D�[?i����9��?�Unknown
�mHostSum"7gradient_tape/binary_crossentropy/logistic_loss/mul/Sum(1�/�$�?9�/�$�?A�/�$�?I�/�$�?ad�D�[?i������?�Unknown
�nHostSum"9gradient_tape/binary_crossentropy/logistic_loss/sub/Sum_1(1�/�$�?9�/�$�?A�/�$�?I�/�$�?ad�D�[?i�RB�,��?�Unknown
woHostReadVariableOp"div_no_nan/ReadVariableOp_1(1X9��v�?9X9��v�?AX9��v�?IX9��v�?a5�]�@?i;�m����?�Unknown
upHostReadVariableOp"div_no_nan/ReadVariableOp(1��Q��?9��Q��?A��Q��?I��Q��?a>w�$�	?i      �?�Unknown*�j
lHostConv2D"sequential/conv1d/conv1d(1B`��"�~@9B`��"�~@AB`��"�~@IB`��"�~@a�ܣk*�?i�ܣk*�?�Unknown
~HostReluGrad"(gradient_tape/sequential/conv1d/ReluGrad(1�O��nVm@9�O��nVm@A�O��nVm@I�O��nVm@a)����?i2����?�Unknown
�HostMaxPoolGrad":gradient_tape/sequential/max_pooling1d/MaxPool/MaxPoolGrad(1����xMm@9����xMm@A����xMm@I����xMm@a��=4���?iV���r�?�Unknown
yHostMatMul"%gradient_tape/sequential/dense/MatMul(1�t��k@9�t��k@A�t��k@I�t��k@a�3�c.m�?i�W��F�?�Unknown
�HostConv2DBackpropFilter";gradient_tape/sequential/conv1d/conv1d/Conv2DBackpropFilter(1�z�G=f@9�z�G=f@A�z�G=f@I�z�G=f@a�S�P?i�sv�c��?�Unknown
uHostMaxPool" sequential/max_pooling1d/MaxPool(1���Q�d@9���Q�d@A���Q�d@I���Q�d@aQ��?iw#����?�Unknown
nHostBiasAdd"sequential/conv1d/BiasAdd(1y�&1�]@9y�&1�]@Ay�&1�]@Iy�&1�]@a�5\��?i�s\�ޙ�?�Unknown
{HostMatMul"'gradient_tape/sequential/dense/MatMul_1(15^�I�T@95^�I�T@A5^�I�T@I5^�I�T@a�g��@Z�?iPzF����?�Unknown
�	HostResourceApplyAdam"$Adam/Adam/update_2/ResourceApplyAdam(1��|?5�P@9��|?5�P@A��|?5�P@I��|?5�P@a23}���?i�k/�Ǭ�?�Unknown
o
Host_FusedMatMul"sequential/dense/Relu(1�l���N@9�l���N@A�l���N@I�l���N@as9C2�T�?i����m�?�Unknown
�HostDataset"3Iterator::Model::ParallelMap::Zip[1]::ForeverRepeat(1)\����?@9)\����?@A%��C�<@I%��C�<@a���cD�?i���b��?�Unknown
�HostBiasAddGrad"3gradient_tape/sequential/conv1d/BiasAdd/BiasAddGrad(11�Z�;@91�Z�;@A1�Z�;@I1�Z�;@a'����y�?iQ�|nfF�?�Unknown
hHostRelu"sequential/conv1d/Relu(1     �7@9     �7@A     �7@I     �7@aD�,^�?i�i��3��?�Unknown
fHostGreaterEqual"GreaterEqual(1}?5^��1@9}?5^��1@A}?5^��1@I}?5^��1@a;/O&��~?i�>m3��?�Unknown
�HostResourceApplyAdam"$Adam/Adam/update_4/ResourceApplyAdam(1�&1��0@9�&1��0@A�&1��0@I�&1��0@av1��}?iO`i1c�?�Unknown
�HostResourceApplyAdam"$Adam/Adam/update_5/ResourceApplyAdam(1���S�e/@9���S�e/@A���S�e/@I���S�e/@aq+ȩ\�{?i���}J�?�Unknown
dHostDataset"Iterator::Model(1�z�G�:@9�z�G�:@A��S�[/@I��S�[/@a2��O`�{?i`�\����?�Unknown
lHostIteratorGetNext"IteratorGetNext(1=
ףp}*@9=
ףp}*@A=
ףp}*@I=
ףp}*@a�*���>w?i��J��?�Unknown
�HostResourceApplyAdam""Adam/Adam/update/ResourceApplyAdam(1V-���(@9V-���(@AV-���(@IV-���(@a������u?ix��eW��?�Unknown
[HostAddV2"Adam/add(1��~j�4(@9��~j�4(@A��~j�4(@I��~j�4(@a�6 ��=u?i����?�Unknown
gHostStridedSlice"strided_slice(1P��n�&@9P��n�&@AP��n�&@IP��n�&@a���V�s?iBV���-�?�Unknown
qHostDataset"Iterator::Model::ParallelMap(1���x�f&@9���x�f&@A���x�f&@I���x�f&@a��y�w�s?iJvU�?�Unknown
^HostGatherV2"GatherV2(1������%@9������%@A������%@I������%@a>�me;!s?iF%��X{�?�Unknown
�HostDataset"=Iterator::Model::ParallelMap::Zip[0]::FlatMap[0]::Concatenate(1�O��n,@9�O��n,@A�V�#@I�V�#@a\ٌ�Dq?i�>Z��?�Unknown
tHost_FusedMatMul"sequential/dense_1/BiasAdd(1�"��~j"@9�"��~j"@A�"��~j"@I�"��~j"@a����	)p?ijl�-3��?�Unknown
tHostAssignAddVariableOp"AssignAddVariableOp(1J+� @9J+�@AJ+� @IJ+�@aV�g~"<l?iG�7Po��?�Unknown
�HostResourceApplyAdam"$Adam/Adam/update_3/ResourceApplyAdam(1H�z��@9H�z��@AH�z��@IH�z��@a�V\-�k?i�0<}R��?�Unknown
oHostReadVariableOp"Adam/ReadVariableOp(1�(\��u@9�(\��u@A�(\��u@I�(\��u@a��!5J�k?i�Rq���?�Unknown
{HostMatMul"'gradient_tape/sequential/dense_1/MatMul(1V-�@9V-�@AV-�@IV-�@a�b�j?im�)��+�?�Unknown
`HostGatherV2"
GatherV2_1(1;�O���@9;�O���@A;�O���@I;�O���@a��:dci?i�,0`E�?�Unknown
|HostSelect"(binary_crossentropy/logistic_loss/Select(1#��~j�@9#��~j�@A#��~j�@I#��~j�@a'��J\7i?i��w��^�?�Unknown
v HostDataset"!Iterator::Model::ParallelMap::Zip(1�Zd;�K@9�Zd;�K@AZd;�O@IZd;�O@a�5Á�g?i�f��Tv�?�Unknown
�!HostReadVariableOp"'sequential/dense/BiasAdd/ReadVariableOp(1���x�@9���x�@A���x�@I���x�@a�'��3cg?i� �Ϸ��?�Unknown
�"HostBiasAddGrad"4gradient_tape/sequential/dense_1/BiasAdd/BiasAddGrad(1�I+�@9�I+�@A�I+�@I�I+�@a�p �ff?i����?�Unknown
v#HostMul"%binary_crossentropy/logistic_loss/mul(1����M�@9����M�@A����M�@I����M�@a��x�e?i�������?�Unknown
�$HostReadVariableOp"4sequential/conv1d/conv1d/ExpandDims_1/ReadVariableOp(1�(\���@9�(\���@A�(\���@I�(\���@a)t�db?i��Y��?�Unknown
Y%HostPow"Adam/Pow(1-���F@9-���F@A-���F@I-���F@ao��a?iA��#��?�Unknown
~&HostSelect"*binary_crossentropy/logistic_loss/Select_1(1-���F@9-���F@A-���F@I-���F@a�ʃ<;�`?i��/��?�Unknown
�'HostResourceApplyAdam"$Adam/Adam/update_1/ResourceApplyAdam(1V-2@9V-2@AV-2@IV-2@ah.��B�`?i:V`r���?�Unknown
�(HostDynamicStitch"/gradient_tape/binary_crossentropy/DynamicStitch(1��K7��@9��K7��@A��K7��@I��K7��@a��īju`?i�[�?�Unknown
V)HostMean"Mean(1���K�@9���K�@A���K�@I���K�@a����_?i�a����?�Unknown
�*HostTile"6gradient_tape/binary_crossentropy/weighted_loss/Tile_1(1j�t��@9j�t��@Aj�t��@Ij�t��@a�!���^?i��{X/�?�Unknown
\+HostGreater"Greater(1�Zd;�@9�Zd;�@A�Zd;�@I�Zd;�@a�؟B`�]?i�B4�&>�?�Unknown
�,HostDataset"MIterator::Model::ParallelMap::Zip[0]::FlatMap[0]::Concatenate[0]::TensorSlice(1�n���@9�n���@A�n���@I�n���@a���nx]?iġ��L�?�Unknown
V-HostSum"Sum_2(1'1�Z@9'1�Z@A'1�Z@I'1�Z@a/�!!��\?iU2B<[�?�Unknown
~.HostMaximum")gradient_tape/binary_crossentropy/Maximum(1V-�@9V-�@AV-�@IV-�@a9����Z?i�Q��h�?�Unknown
�/HostSelect"8gradient_tape/binary_crossentropy/logistic_loss/Select_3(1����K@9����K@A����K@I����K@a���͕Z?i����u�?�Unknown
�0Host	ZerosLike":gradient_tape/binary_crossentropy/logistic_loss/zeros_like(1�A`��"@9�A`��"@A�A`��"@I�A`��"@a�FS��qZ?i�s;�7��?�Unknown
z1HostLog1p"'binary_crossentropy/logistic_loss/Log1p(1`��"��@9`��"��@A`��"��@I`��"��@a��-�MZ?i�N��^��?�Unknown
}2HostMatMul")gradient_tape/sequential/dense_1/MatMul_1(1�V-@9�V-@A�V-@I�V-@a��'4�Y?i� ��+��?�Unknown
]3HostCast"Adam/Cast_1(1���K7@9���K7@A���K7@I���K7@a�����X?i�A���?�Unknown
j4HostMean"binary_crossentropy/Mean(1�rh��|	@9�rh��|	@A�rh��|	@I�rh��|	@a�Ʒ]V?if_����?�Unknown
V5HostAddN"AddN(1�� �rh	@9�� �rh	@A�� �rh	@I�� �rh	@a&��KV?iy�/���?�Unknown
�6HostDataset"?Iterator::Model::ParallelMap::Zip[1]::ForeverRepeat::FromTensor(1!�rh��@9!�rh���?A!�rh��@I!�rh���?a{����U?i�������?�Unknown
v7HostExp"%binary_crossentropy/logistic_loss/Exp(1H�z�G@9H�z�G@AH�z�G@IH�z�G@a��r+YNU?i�f3y��?�Unknown
j8HostReadVariableOp"ReadVariableOp(1�Q���@9�Q���@A�Q���@I�Q���@a##�:�T?i�Љ���?�Unknown
�9HostReadVariableOp")sequential/dense_1/BiasAdd/ReadVariableOp(1�G�z�@9�G�z�@A�G�z�@I�G�z�@a�MN��S?i4k`�z��?�Unknown
V:HostCast"Cast(1��Q��@9��Q��@A��Q��@I��Q��@a�����S?i�C<|\��?�Unknown
X;HostEqual"Equal(1��Q��@9��Q��@A��Q��@I��Q��@a�����S?iH�=��?�Unknown
�<HostDataset"-Iterator::Model::ParallelMap::Zip[0]::FlatMap(1V-�0@9V-�0@AH�z�G@IH�z�G@a2?E��S?i�>��?�Unknown
�=HostBiasAddGrad"2gradient_tape/sequential/dense/BiasAdd/BiasAddGrad(1�z�G�@9�z�G�@A�z�G�@I�z�G�@a�0<43S?i ����?�Unknown
X>HostCast"Cast_2(1��(\��@9��(\��@A��(\��@I��(\��@a֡�C�7R?i����?�Unknown
�?HostAddV2"3gradient_tape/binary_crossentropy/logistic_loss/add(1q=
ףp@9q=
ףp@Aq=
ףp@Iq=
ףp@a�/�t��Q?ii/�ñ"�?�Unknown
�@Host	ZerosLike"<gradient_tape/binary_crossentropy/logistic_loss/zeros_like_1(1q=
ףp@9q=
ףp@Aq=
ףp@Iq=
ףp@a�/�t��Q?i����+�?�Unknown
�AHostGreaterEqual".binary_crossentropy/logistic_loss/GreaterEqual(1��Q�@9��Q�@A��Q�@I��Q�@a�z�ԧQ?i`G܉}4�?�Unknown
oBHostSigmoid"sequential/dense_1/Sigmoid(1�n���@9�n���@A�n���@I�n���@aG<E�%^Q?i�i��,=�?�Unknown
vCHostNeg"%binary_crossentropy/logistic_loss/Neg(1��ʡE�@9��ʡE�@A��ʡE�@I��ʡE�@aȟv-LQ?iN�F��E�?�Unknown
vDHostAssignAddVariableOp"AssignAddVariableOp_2(1��n��@9��n��@A��n��@I��n��@a�IǢ�tP?i���N�?�Unknown
�EHostBroadcastTo"-gradient_tape/binary_crossentropy/BroadcastTo(1�(\���@9�(\���@A�(\���@I�(\���@a1y`]3�O?iaoB�U�?�Unknown
~FHostAssignAddVariableOp"Adam/Adam/AssignAddVariableOp(1o��ʡ@9o��ʡ@Ao��ʡ@Io��ʡ@a�u����N?inַ�]�?�Unknown
�GHostFloorDiv"*gradient_tape/binary_crossentropy/floordiv(1F����x@9F����x@AF����x@IF����x@a����N?i/�4Ue�?�Unknown
}HHostReluGrad"'gradient_tape/sequential/dense/ReluGrad(1�O��n@9�O��n@A�O��n@I�O��n@a��;�M?iiwz��l�?�Unknown
�IHostSelect"8gradient_tape/binary_crossentropy/logistic_loss/Select_2(1P��n� @9P��n� @AP��n� @IP��n� @a(X\��L?i��kt�?�Unknown
tJHostReadVariableOp"Adam/Cast/ReadVariableOp(1��C�l��?9��C�l��?A��C�l��?I��C�l��?aB��K?iqR�-{�?�Unknown
xKHostCast"&gradient_tape/binary_crossentropy/Cast(1�$��C�?9�$��C�?A�$��C�?I�$��C�?aQ�|BoK?i*y����?�Unknown
�LHost
Reciprocal":gradient_tape/binary_crossentropy/logistic_loss/Reciprocal(1��/�$�?9��/�$�?A��/�$�?I��/�$�?akVN��sJ?i�L�艈�?�Unknown
�MHostSelect"6gradient_tape/binary_crossentropy/logistic_loss/Select(1m������?9m������?Am������?Im������?at���+J?i�Q<���?�Unknown
bNHostDivNoNan"div_no_nan_1(1+���?9+���?A+���?I+���?a=onF�PI?iU��h��?�Unknown
wOHostReadVariableOp"div_no_nan_1/ReadVariableOp(1J+��?9J+��?AJ+��?IJ+��?ag���%�G?i��EZ��?�Unknown
vPHostSub"%binary_crossentropy/logistic_loss/sub(1���S��?9���S��?A���S��?I���S��?a��%.��F?i"N����?�Unknown
�QHostSum"5gradient_tape/binary_crossentropy/logistic_loss/Sum_1(1����S�?9����S�?A����S�?I����S�?a���^�9F?i��� ���?�Unknown
�RHostMul"3gradient_tape/binary_crossentropy/logistic_loss/mul(1㥛� ��?9㥛� ��?A㥛� ��?I㥛� ��?a��v� �E?i�����?�Unknown
XSHostCast"Cast_3(1��ʡE�?9��ʡE�?A��ʡE�?I��ʡE�?a�/J�C�C?i*��qְ�?�Unknown
rTHostAdd"!binary_crossentropy/logistic_loss(1��ʡE�?9��ʡE�?A��ʡE�?I��ʡE�?a�/J�C�C?i��B���?�Unknown
[UHostPow"
Adam/Pow_1(1+�����?9+�����?A+�����?I+�����?a����FC?i-����?�Unknown
XVHostCast"Cast_4(1�ʡE���?9�ʡE���?A�ʡE���?I�ʡE���?a���`CC?i�A�[��?�Unknown
vWHostSum"%binary_crossentropy/weighted_loss/Sum(1� �rh��?9� �rh��?A� �rh��?I� �rh��?a����+A?i��4Φ��?�Unknown
�XHostReadVariableOp"(sequential/dense_1/MatMul/ReadVariableOp(1Zd;�O��?9Zd;�O��?AZd;�O��?IZd;�O��?a�f٨;(A?iX����?�Unknown
�YHostReadVariableOp"(sequential/conv1d/BiasAdd/ReadVariableOp(1sh��|?�?9sh��|?�?Ash��|?�?Ish��|?�?a ����@?i�t^�)��?�Unknown
�ZHostReadVariableOp"&sequential/dense/MatMul/ReadVariableOp(1��� �r�?9��� �r�?A��� �r�?I��� �r�?a3���90@?i[U��5��?�Unknown
v[HostReadVariableOp"Adam/Cast_2/ReadVariableOp(1;�O��n�?9;�O��n�?A;�O��n�?I;�O��n�?a�׌ӡ,@?i�8QA��?�Unknown
~\HostRealDiv")gradient_tape/binary_crossentropy/truediv(1Zd;�O�?9Zd;�O�?AZd;�O�?IZd;�O�?a���b>?i��PR��?�Unknown
v]HostReadVariableOp"Adam/Cast_3/ReadVariableOp(1�x�&1�?9�x�&1�?A�x�&1�?I�x�&1�?a>���"<?iC�����?�Unknown
�^HostCast"3binary_crossentropy/weighted_loss/num_elements/Cast(1m������?9m������?Am������?Im������?at���+:?i� <*���?�Unknown
�_HostMul"5gradient_tape/binary_crossentropy/logistic_loss/mul_1(1w��/��?9w��/��?Aw��/��?Iw��/��?a��djT9?iRm���?�Unknown
v`HostAssignAddVariableOp"AssignAddVariableOp_3(1�$��C�?9�$��C�?A�$��C�?I�$��C�?a�T@^��7?i]5�C���?�Unknown
�aHostNeg"7gradient_tape/binary_crossentropy/logistic_loss/sub/Neg(1�x�&1�?9�x�&1�?A�x�&1�?I�x�&1�?ab���5?i��4���?�Unknown
TbHostMul"Mul(1bX9���?9bX9���?AbX9���?IbX9���?a	H�3?i�6�!��?�Unknown
�cHostMul"7gradient_tape/binary_crossentropy/logistic_loss/mul/Mul(1�n����?9�n����?A�n����?I�n����?a��rSo3?i� υ��?�Unknown
vdHostAssignAddVariableOp"AssignAddVariableOp_1(1�E�����?9�E�����?A�E�����?I�E�����?a�.���G2?iK?�����?�Unknown
�eHostSum"3gradient_tape/binary_crossentropy/logistic_loss/Sum(1P��n��?9P��n��?AP��n��?IP��n��?a�����1?ic����?�Unknown
`fHostDivNoNan"
div_no_nan(1�Zd;��?9�Zd;��?A�Zd;��?I�Zd;��?a��xp1?i�O�<��?�Unknown
�gHostDivNoNan"@gradient_tape/binary_crossentropy/weighted_loss/value/div_no_nan(1bX9���?9bX9���?AbX9���?IbX9���?a\:�En+?io3����?�Unknown
�hHostSum"7gradient_tape/binary_crossentropy/logistic_loss/sub/Sum(1w��/��?9w��/��?Aw��/��?Iw��/��?a��djT)?i[
A���?�Unknown
�iHostNeg"3gradient_tape/binary_crossentropy/logistic_loss/Neg(1��Q���?9��Q���?A��Q���?I��Q���?avI��1&(?i��#���?�Unknown
}jHostDivNoNan"'binary_crossentropy/weighted_loss/value(1�z�G��?9�z�G��?A�z�G��?I�z�G��?a�e-�k�'?iV��
~��?�Unknown
ykHostReadVariableOp"div_no_nan_1/ReadVariableOp_1(1�/�$�?9�/�$�?A�/�$�?I�/�$�?a�6��{�%?i��b���?�Unknown
�lHostSum"7gradient_tape/binary_crossentropy/logistic_loss/mul/Sum(1�/�$�?9�/�$�?A�/�$�?I�/�$�?a�6��{�%?i<WZ�<��?�Unknown
�mHostSum"9gradient_tape/binary_crossentropy/logistic_loss/sub/Sum_1(1�/�$�?9�/�$�?A�/�$�?I�/�$�?a�6��{�%?i�����?�Unknown
wnHostReadVariableOp"div_no_nan/ReadVariableOp_1(1X9��v�?9X9��v�?AX9��v�?IX9��v�?a!��9e�#?i�.ix���?�Unknown
uoHostReadVariableOp"div_no_nan/ReadVariableOp(1��Q��?9��Q��?A��Q��?I��Q��?abmy�"?i�������?�Unknown