import csv

def convert():
    input_file = 'results/differential_all_variants.csv'
    output_file = 'input_for_vep.vcf'
    
    try:
        with open(input_file, 'r', encoding='utf-8') as f_in, \
             open(output_file, 'w', encoding='utf-8') as f_out:
            
            reader = csv.DictReader(f_in)
            
            # 1. VCF 필수 헤더 작성 (VEP 인식용)
            f_out.write("##fileformat=VCFv4.2\n")
            f_out.write("##source=DifferentialVariantAnalysis\n")
            f_out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            
            count = 0
            for row in reader:
                # CSV 컬럼명에 맞춰 데이터 추출
                chrom = row['chrom']
                pos = row['pos']
                ref = row['ref']
                alt = row['alt']
                
                # INFO 필드에 나중에 참고할 재발률(Recurrence) 정보만 살짝 넣어줍니다
                info = f"RECUR={row['recurrence']};DAR={row['delta_alt_ratio']}"
                
                # VCF 한 줄 작성 (ID, QUAL, FILTER는 '.' 또는 'PASS'로 채움)
                f_out.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t{info}\n")
                count += 1
                
        print(f"✅ 변환 완료! {output_file} 파일이 생성되었습니다. (총 {count}개 변이)")
        
    except FileNotFoundError:
        print(f"❌ 에러: {input_file} 파일을 찾을 수 없습니다. 경로를 확인해주세요.")

if __name__ == "__main__":
    convert()